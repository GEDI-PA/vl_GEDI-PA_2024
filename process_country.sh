#!/bin/bash

# Ensure script exits if any command fails
set -e

# Function to display help message
function show_help {
    echo "Usage: $0 -c <ISO3_country_code> -w <GEDI_week>"
    echo ""
    echo "Options:"
    echo "  -c  ISO3 country code (e.g., ECU for Ecuador)"
    echo "  -w  GEDI week (e.g., 24)"
    echo "  -h  Show this help message"
}

# Parse command-line arguments
while getopts ":c:w:h" opt; do
  case ${opt} in
    c )
      iso3=$OPTARG
      ;;
    w )
      gediwk=$OPTARG
      ;;
    h )
      show_help
      exit 0
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      show_help
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      show_help
      exit 1
      ;;
  esac
done

# Check if ISO3 country code and GEDI week are set
if [ -z "${iso3}" ] || [ -z "${gediwk}" ]; then
    echo "Both ISO3 country code and GEDI week are required."
    show_help
    exit 1
fi

# Set the R script path
R_SCRIPT="global_process_part1_2024_MAAP_step123.R"

# Execute the R script with the provided ISO3 country code and GEDI week
Rscript $R_SCRIPT $iso3 $gediwk
