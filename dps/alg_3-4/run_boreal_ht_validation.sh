#!/bin/bash

libdir=$(dirname "$(dirname "${basedir}")")/lib

mkdir output
source activate icesat2_boreal

args=(--boreal_tiles_path "${1}")
args+=(--ht_tindex_path "${2}")
args+=(--lvis_footprints_path "${3}")
args+=(--boreal_tile_num "${4}")
args+=(--metrics "${5}")
args+=(--year "${6}")

command=(python "${libdir}/ht_validation.py" "${args[@]}")
echo "${command[@]}"
"${command[@]}"
