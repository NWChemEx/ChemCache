#!/bin/sh

# This script downloads reference data.
#
# Usage:
#   download_reference_data.sh
#
# Arguments:
#   doxygen_target_name: The name of the CMake target which controls building
#                        the Doxygen documentation. The name of this
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

# Create necessary data subdirectories
mkdir -p "${BASIS_SETS}"
#mkdir -p "${DENSITIES}"
#mkdir -p "${MOLECULES}"

# Activate virtual environment
. venv/bin/activate

# Call download script(s)

echo "Calling ${REFERENCE_DATA}/scrape_bse.py ${BASIS_SETS}"
${PYTHON} ${REFERENCE_DATA}/scrape_bse.py ${BASIS_SETS}
