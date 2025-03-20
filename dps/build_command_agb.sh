#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
source activate r
conda install -c conda-forge r-optparse -y
