#!/bin/bash
# this is intended for running DPS jobs; the input directory is where four files have been pulled because download=TRUE in the algorithm_config.yaml file
# a tar file of biomass models, a data table csv, and two raster stack geotiff files

#source activate icesat2_boreal
basedir=$( cd "$(dirname "$0")" ; pwd -P )
libdir=$(dirname "$(dirname "${basedir}")")/lib
#unset PROJ_LIB

#pip install --user -r ${basedir}/requirements.txt

mkdir output

# Note: the numbered args are fed in with the in_param_dict in the Run DPS chunk of 3.4_dps.ipynb
ATL08_tindex_master_fn=${1}
ATL08_SAMPLE_CSV=${2}
in_tile_num=${3}
in_tile_fn=${4}
in_tile_field=${5}
minDOY=${6}
maxDOY=${7}
max_sol_el=${8}
expand_training=${9}
local_train_perc=${10}
min_n=${11}
predict_var=${12}
max_n=${13}
pred_vars=${14}
bio_models_tar_fn=${15}
year=${16}

source activate icesat2_boreal

tar -xf ${bio_models_tar_fn}

# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

# Get the output merged CSV of filtered ATL08 for the input tile and its neighbors
cmd="python ${libdir}/merge_neighbors_atl08.py -in_tile_num ${in_tile_num} -in_tile_fn ${in_tile_fn} -in_tile_field ${in_tile_field} -csv_list_fn ${ATL08_tindex_master_fn} -out_dir ${OUTPUTDIR}"

echo $cmd
eval $cmd

# Set the output merged CSV name to a var
MERGED_ATL08_CSV=$(ls ${OUTPUTDIR}/*_merge_neighbors_*.csv | head -1)

echo $MERGED_ATL08_CSV
echo $ATL08_SAMPLE_CSV

source activate r

# required arguments
args=(--atl08_path "${MERGED_ATL08_CSV}")
args+=(--broad_path "${2}")
args+=(--year "${16}")

# optional arguments
[[ -n "${6}" ]] && args+=(--minDOY "${6}")
[[ -n "${7}" ]] && args+=(--maxDOY "${7}")
[[ -n "${8}" ]] && args+=(--max_sol_el "${8}")
[[ -n "${9}" ]] && args+=(--expand_training "${9}")
[[ -n "${10}" ]] && args+=(--local_train_perc "${10}")
[[ -n "${11}" ]] && args+=(--min_samples "${11}")
[[ -n "${13}" ]] && args+=(--max_samples "${13}")
[[ -n "${12}" ]] && args+=(--predict_var "${12}")
[[ -n "${14}" ]] && args+=(--pred_vars "${14}")

command=(Rscript "${libdir}/model_comparison.R" "${args[@]}")
echo "${command[@]}"
"${command[@]}"
