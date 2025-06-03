#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

#install requirements packages
conda env create -f ${basedir}/above_env_r_3.1.4.yml

pushd ${HOME}

source activate r
conda install -c conda-forge r-optparse r-ranger r-ggplot2 -y
