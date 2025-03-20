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
LIDAR_CSV=${1}
TOPO_TIF=${2}
HLS_TIF=${3}
LC_TIF=${4}
DO_SLOPE_VALID_MASK=${5}
ATL08_SAMPLE_CSV=${6}
iters=${7}
calculate_uncertainty=${8}
minDOY=${9}
maxDOY=${10}
max_sol_el=${11}
expand_training=${12}
local_train_perc=${13}
min_n=${14}
boreal_vect_fn=${15}
predict_var=${16}
max_n=${17}
pred_vars=${18}
bio_models_tar_fn=${19}
SAR_TIF=${20}
year=${21}

tar -xf ${bio_models_tar_fn}

# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

source activate r

# required arguments
# this could actaully be lidar data in csv format from atl08, GEDI, or mix of both
args=(--atl08_path "${1}")
args+=(--broad_path "${6}")
args+=(--topo_path "${2}")
args+=(--hls_path "${3}")
args+=(--lc_path "${4}")
args+=(--boreal_vector_path "${15}")
args+=(--year "${21}")

# optional arguments
[[ -n "${20}" ]] && args+=(--sar_path "${20}")
[[ -n "${5}" ]] && args+=(--mask "${5}")
[[ -n "${7}" ]] && args+=(--uncertainty_iterations "${7}")
[[ -n "${8}" ]] && args+=(--calculate_uncertainty "${8}")
[[ -n "${9}" ]] && args+=(--minDOY "${9}")
[[ -n "${10}" ]] && args+=(--maxDOY "${10}")
[[ -n "${11}" ]] && args+=(--max_sol_el "${11}")
[[ -n "${12}" ]] && args+=(--expand_training "${12}")
[[ -n "${13}" ]] && args+=(--local_train_perc "${13}")
[[ -n "${14}" ]] && args+=(--min_samples "${14}")
[[ -n "${17}" ]] && args+=(--max_samples "${17}")
[[ -n "${16}" ]] && args+=(--predict_var "${16}")
[[ -n "${18}" ]] && args+=(--pred_vars "${18}")

command=(Rscript "${libdir}/mapBoreal_simple.R" "${args[@]}")
echo "${command[@]}"
"${command[@]}"
