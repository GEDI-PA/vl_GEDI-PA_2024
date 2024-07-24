#!/bin/bash


while getopts ":c:w:h" opt; do
  case ${opt} in
    c )
      iso3=$OPTARG
      ;;
    w )
      gediwk=$OPTARG
      ;;
  esac
done

# Set the R script path
R_SCRIPT="global_process_part1_2024_MAAP_step123.R"

# Execute the R script with the provided ISO3 country code and GEDI week
Rscript global_process_part1_2024_MAAP_step123.R $iso3 $gediwk
