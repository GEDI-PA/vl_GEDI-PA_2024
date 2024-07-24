#!/bin/bash

# Initialize variables
country=""
wk=""

# Parse command-line options
while getopts c:w: flag
do
	case "${flag}" in
		c) country=${OPTARG};;
		w) wk=${OPTARG};;
	esac
done

# Check if necessary options are provided
if [[ -z "$country" || -z "$wk" ]]; then
	echo "Usage: $0 -c <country> -w <week>"
	exit 1
fi

# Run the R script for the specified country and week
Rscript global_process_part1_2024_MAAP_step123.R $country $wk
