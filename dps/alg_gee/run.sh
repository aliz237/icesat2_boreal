#!/bin/bash
# this is intended for running DPS jobs; the input directory is where a single file has been pulled because download=TRUE in the algorithm_config.yaml file

# This installs the python libs needed to run the script at the bottom
# these libs are NOT included in the base image (vanilla: https://mas.maap-project.org/root/ade-base-images/-/blob/vanilla/docker/Dockerfile)
#conda install -yq -c conda-forge geopandas rio-cogeo rio-tiler importlib_resources
set -x
source activate icesat2_boreal

unset PROJ_LIB
#conda list
#pip install --user -U numpy==1.20.3 geopandas==0.9.0 rio-cogeo==2.3.1 rio-tiler==2.1.4 rasterio==1.2.6 morecantile==2.1.4 pystac-client importlib_resources 

mkdir output

basedir=$( cd "$(dirname "$0")" ; pwd -P )  # goes to alg_3-1-5/

# First file in input/ dir
# TODO: Fragile relying on alphabetical order

# DPS notebook used to submit this job will make a new credentials file and place it in the private bucket on s3 at the 
INPUT1='credentials'

mkdir -p .config/earthengine/
cp ${INPUT1} .config/earthengine/

## Hard coded args for each run (if any; usually just output dir)

# Work dir is always from where your script is called
# Base dir is always the relative dir within the run*.sh script

# Absolute path here
# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

python ${basedir}/../../lib/export_gee_to_maap.py \
--in_tile_num ${1} \
--dims ${2} \
--asset_path ${3} \
--out_dir ${OUTPUTDIR}