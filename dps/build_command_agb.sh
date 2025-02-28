#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
echo ${basedir}
#install requirements packages
conda env create -f ${basedir}/above_env_r_3.1.4.yml

pushd ${HOME}

source activate r
conda env update -f ${basedir}/r_libs.yml
Rscript -e "install.packages('paws.storage', repos='https://cloud.r-project.org')"
