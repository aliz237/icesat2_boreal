#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

conda env create -f ${basedir}/above_env_r_3.1.4.yml
pushd ${HOME}

source activate icesat2_boreal
