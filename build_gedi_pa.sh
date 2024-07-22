#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

#install requirements packages
conda env update -f ${basedir}/gedi_pa_env.yml

pushd ${HOME}

source activate gedi_pa_env

