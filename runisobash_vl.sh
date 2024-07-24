#!/bin/bash
basedir=$( cd "$(dirname "$0")" ; pwd -P )

# Initialize variables
country=""

# Parse command-line options
while getopts c: flag
do
	case "${flag}" in
		c) country=${OPTARG};;
	esac
done

Rscript ${basedir}/global_process_part1_2024_MAAP_step123.R $country