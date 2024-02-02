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

# Set reference data directory
REF_DATA="reference_data"

# Set data directories
ATOMIC_INFO="${REF_DATA}/physical_data"
BASIS_SETS="${REF_DATA}/basis_sets/default"
DENSITIES="${REF_DATA}/atomic_densities/default"
MOLECULES="${REF_DATA}/molecules/default"

# Activate virtual environment
. venv/bin/activate

# Call download script(s)
SCRIPTS_DIR="utils/data_management"

echo "Calling ${SCRIPTS_DIR}/scrape_bse.py ${BASIS_SETS}"
python ${SCRIPTS_DIR}/scrape_bse.py ${BASIS_SETS}
