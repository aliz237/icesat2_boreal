from collections import defaultdict
import argparse

import numpy as np
import geopandas as gpd
import pandas as pd

import rasterio
from shapely.geometry import box, Point

from osgeo import gdal

import s3fs
import richdem as rd


def get_datasets(tile_num, ht_tindex_path, lvis_footprints_path, boreal_tiles_path, metrics):
    ht_tindex = pd.read_csv(ht_tindex_path)
    boreal_ht_raster_path = ht_tindex.set_index('tile_num').loc[tile_num,'s3_path']

    lvis_tiles = gpd.read_file(lvis_footprints_path)

    boreal_tiles = gpd.read_file(boreal_tiles_path)
    boreal_tiles_reproj = boreal_tiles.to_crs(lvis_tiles.crs)

    bbox_poly = box(*lvis_tiles.total_bounds)
    boreal_tiles_reproj_filt = boreal_tiles_reproj[boreal_tiles_reproj.intersects(bbox_poly)]

    lvis_paths = download_intersecting_lvis_tiles(
        tile_num,
        boreal_tiles_reproj_filt.set_index('tile_num'),
        lvis_tiles,
        metrics
    )
    return boreal_ht_raster_path, lvis_paths


def download_intersecting_lvis_tiles(boreal_tile_num, boreal_tiles, lvis_tiles, metrics):
    boreal_tile_geom = boreal_tiles.loc[boreal_tile_num, 'geometry']
    lvis_cnt_paths = lvis_tiles[lvis_tiles.intersects(boreal_tile_geom)][['FILE']].values.flatten()
    s3_base = 's3://ornl-cumulus-prod-protected/above/ABoVE_LVIS_VegetationStructure/data/'
    lvis_s3_paths = {
        metric:
        [s3_base + y.replace('lvis_pt_cnt_30m.tif', '') + metric + '_30m.tif' for y in lvis_cnt_paths]
        for metric in metrics
    }
    fs = s3fs.S3FileSystem()
    lvis_local_paths = defaultdict(list)
    for metric, s3_paths in lvis_s3_paths.items():
        for s3_path in s3_paths:
            local_path = s3_path.replace(s3_base, '/tmp/')
            lvis_local_paths[metric].append(local_path)
            fs.get(s3_path, local_path)
            print(f"Downloaded {s3_path} to {local_path}")
    return lvis_local_paths


def raster_to_gdf(raster_path, value_column='value'):
    with rasterio.open(raster_path) as src:
        data = src.read(1)
        transform = src.transform
        nodata = src.nodata
        crs = src.crs
        print(raster_path, nodata, crs, src.shape)

        rows, cols = np.where(data != nodata)
        xs, ys = rasterio.transform.xy(transform, rows, cols)

        gdf = gpd.GeoDataFrame(
            {value_column: data[rows, cols]},
            geometry=[Point(x, y) for x, y in zip(xs, ys)],
            crs=crs
        )

    return gdf


def lvis_boreal_vrt(base_meta, lvis_raster_paths, metric):
    print('Building the VRT...')
    vrt = gdal.BuildVRT("/vsimem/mosaic.vrt", lvis_raster_paths)

    print('Warping the VRT to match the boreal raster...')
    warp_options = gdal.WarpOptions(
        format="GTiff",
        outputBounds=base_meta['bounds'],
        width=base_meta['width'],
        height=base_meta['height'],
        dstSRS=base_meta['crs'].to_wkt(),
        resampleAlg="near"
    )
    out_vrt_fn = f"/vsimem/{metric}.vrt"
    out = gdal.Warp(out_vrt_fn, vrt, options=warp_options)
    if out is None:
        raise RuntimeError("Warp operation failed.")

    return out_vrt_fn


def sample_raster_at_points(gdf, raster_path, column_name):
    with rasterio.open(raster_path) as src:
        coords = [(geom.x, geom.y) for geom in gdf.geometry]
        sampled_values = list(src.sample(coords))
        sampled_values = [val[0] for val in sampled_values]
        gdf[column_name] = sampled_values
    return gdf


def rdarray_to_vsimem(rdarr, crs):
    height, width = rdarr.shape
    if height == 0 or width == 0:
        raise ValueError("Input rdarray has invalid shape (0 dimensions).")
    print(rdarr.shape)

    if not rdarr.geotransform or len(rdarr.geotransform) != 6:
        raise ValueError("Invalid geotransform provided.")
    print(rdarr.geotransform)

    driver = gdal.GetDriverByName("GTiff")
    vsimem_path = "/vsimem/slope_raster.tif"

    dst_ds = driver.Create(vsimem_path, width, height, 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform(rdarr.geotransform)
    dst_ds.SetProjection(crs)
    dst_ds.GetRasterBand(1).WriteArray(rdarr)
    dst_ds.GetRasterBand(1).SetNoDataValue(rdarr.no_data)
    dst_ds.FlushCache()
    dst_ds = None

    return vsimem_path


def base_metadata(boreal_ht_raster_path):
    with rasterio.open(boreal_ht_raster_path) as base_raster:
        base_meta = base_raster.meta
        base_meta['bounds'] = base_raster.bounds
        base_meta['transform'] = base_raster.transform
    return base_meta


def slope_from_ZG_mean(ZG_mean_rsater_path, base_meta):
    lvis_ZG_mean = rd.LoadGDAL(ZG_mean_rsater_path, no_data=255)
    lvis_slope = rd.TerrainAttribute(lvis_ZG_mean, attrib="slope_degrees")
    vsimem_path = rdarray_to_vsimem(lvis_slope, base_meta['crs'].to_wkt())
    return vsimem_path


def join_boreal_ht_to_lvis_as_dataframe(mosaic_dfs, metric_dfs, boreal_ht_raster_path):
    print('joining boreal Ht to first LVIS dataframe...')
    joined = sample_raster_at_points(metric_dfs[0], boreal_ht_raster_path, 'Ht')
    print(joined.head())
    for i, df in enumerate(metric_dfs[1:]):
        print(df.head())
        joined = gpd.sjoin(joined, df, how="inner", predicate="intersects")
        joined = joined.drop(columns=["index_right"])
        print(joined.head())
    return joined


def main(boreal_tile_num, ht_tindex_path, boreal_tiles_path, lvis_footprints_path, metrics, year):
    if "ZG_mean" not in metrics:
        # this is required for slope calculation
        metrics.append('ZG_mean')

    boreal_ht_raster_path, lvis_paths = get_datasets(
        boreal_tile_num, ht_tindex_path, lvis_footprints_path, boreal_tiles_path, metrics
    )

    base_meta = base_metadata(boreal_ht_raster_path)

    mosaic_paths = {
        metric: lvis_boreal_vrt(base_meta, lvis_paths[metric], metric)
        for metric in metrics
    }

    metric_dfs = [raster_to_gdf(mosaic_paths[metric], metric)
                  for metric in metrics]

    joined = join_boreal_ht_to_lvis_as_dataframe(
        mosaic_paths, metric_dfs, boreal_ht_raster_path
    )

    slope_raster_path = slope_from_ZG_mean(mosaic_paths['ZG_mean'], base_meta)

    joined = sample_raster_at_points(joined, slope_raster_path, 'slope')
    joined.to_csv(f'output/Ht_validation_{boreal_tile_num}.csv')


def parse_arguments():
    parser = argparse.ArgumentParser(description="Parser for ht_validation.py arguments.")

    parser.add_argument("--boreal_tiles_path", required=True, help="Path to the boreal tiles file")
    parser.add_argument("--ht_tindex_path", required=True, help="Path to the ht tindex file")
    parser.add_argument("--lvis_footprints_path", required=True, help="Path to the LVIS footprints file")
    parser.add_argument("--boreal_tile_num", required=True, type=int, help="Boreal tile number")
    parser.add_argument("--metrics", required=True, help="Metrics to process")
    parser.add_argument("--year", required=True, type=int, help="Year of the lvis data")

    args = vars(parser.parse_args())
    args['metrics'] = args['metrics'].strip().split()
    print(args)
    print("Parsed arguments:")
    for arg, value in args.items():
        print(f"{arg}: {value}")

    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(**args)
    
