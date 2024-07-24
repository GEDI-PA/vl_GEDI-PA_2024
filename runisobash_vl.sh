#!/bin/bash

# Initialize variables
country=""

# Parse command-line options
while getopts c: flag
do
	case "${flag}" in
		c) country=${OPTARG};;
	esac
done

# Run the R script for the specified country
Rscript global_process_part1_2024_MAAP_step123.R $country
