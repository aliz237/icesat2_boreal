#!/bin/bash
# this is intended for running DPS jobs; the input directory is where a single file has been pulled because download=TRUE in the algorithm_config.yaml file

# This installs the python libs needed to run the script at the bottom
# these libs are NOT included in the base image (vanilla: https://mas.maap-project.org/root/ade-base-images/-/blob/vanilla/docker/Dockerfile)
#conda install -yq -c conda-forge geopandas rio-cogeo rio-tiler importlib_resources
source activate icesat2_boreal
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

unset PROJ_LIB

#install requirements packages
#pip install --user -r ${basedir}/requirements.txt

mkdir output

FILENAME=$(ls -d input/*)

# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

#command line call to try for testing
# python /projects/icesat2_boreal/lib/extract_filter_atl08.py --no-filter-qual --do_30m -i "/projects/test_data/test_data_30m/ATL08_30m_20181014001049_02350102_003_01.h5"

# Call for spring 2022 version
#python ${basedir}/../../lib/extract_filter_atl08.py --i ${FILENAME} --no-filter-qual --do_30m -o output/

# Call for summer/fall 2022 version
#python ${basedir}/../../lib/extract_filter_atl08.py --i ${FILENAME} --no-filter-qual --do_30m --get_gedi_rh -o output/

# Call for summer/fall 2022 version 2
python ${basedir}/../../lib/extract_filter_atl08.py \
--i ${FILENAME} \
--no-filter-qual \
--rh_type ${1} \
--seg_length ${2} \
--out_dir output/
