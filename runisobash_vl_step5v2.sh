#!/bin/bash
basedir=$( cd "$(dirname "$0")" ; pwd -P )

mkdir -p output

# Initialize variables
country=""
flag=""

# Parse command-line options
#while getopts c: flag
#do
#	case "${flag}" in
#		c) country=${OPTARG};;
#	esac
#done
country=${1}
flag=${2}
#output_filename=${2}

#Rscript ${basedir}/global_process_part1_2024_MAAP_step123.R $country
conda run --no-capture-output --name gedi_pa_env Rscript ${basedir}/global_process_part2_2024_MAAP_step5v2.R $country $flag
#--output_file output/${output_filename}