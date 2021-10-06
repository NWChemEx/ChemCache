#!/bin/sh

# This script downloads reference data.
#
# Usage:
#   download_reference_data.sh
#
# Notes:
#   - This script assumes the dependencies in reference_data/requirements.txt 
#     have been installed using PIP into a virtual environment called "venv"

set -e # Exit with error if any command fails

# Set some environment-specific values
PYTHON="python"

# Set reference data directory
REFERENCE_DATA="reference_data"

# Set data directories
ATOMIC_INFO="${REFERENCE_DATA}/physical_data"
BASIS_SETS="${REFERENCE_DATA}/basis_sets/default"
DENSITIES="${REFERENCE_DATA}/atomic_densities/default"
MOLECULES="${REFERENCE_DATA}/molecules/default"

# Remove existing downloaded data files
# rm -rf "${ATOMIC_INFO}"
rm -rf "${BASIS_SETS}"
# rm -rf "${DENSITIES}"
# rm -rf "${MOLECULES}"

# Recreate necessary data subdirectories
# mkdir -p "${ATOMIC_INFO}"
mkdir -p "${BASIS_SETS}"
# mkdir -p "${DENSITIES}"
# mkdir -p "${MOLECULES}"

# Activate virtual environment
. venv/bin/activate

# Call download script(s)

echo "Calling ${REFERENCE_DATA}/scrape_bse.py ${BASIS_SETS}"
${PYTHON} ${REFERENCE_DATA}/scrape_bse.py ${BASIS_SETS}
