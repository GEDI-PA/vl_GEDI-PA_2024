#!/bin/bash
basedir=$( cd "$(dirname "$0")" ; pwd -P )

mkdir -p output
outdir=${PWD}/output

country=${1}

range=${2}


#Rscript ${basedir}/global_process_part1_2024_MAAP_step123.R $country $outdir
cd ${basedir}
conda run --no-capture-output --name gedi_pa_env Rscript global_process_part3_2024_MAAP_step5.r $country $outdir $range 