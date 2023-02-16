#!/bin/sh

# This script generates the reference data source code.
#
# Usage:
#   generate_reference_data.sh
#
# Notes:
#   - This script assumes Sphinx is installed via PIP into a virtual environment
#     called "venv"

set -e # Exit with error if any command fails

# Set some environment-specific values
PYTHON="python"

# Set reference data directory
REFERENCE_DATA="reference_data"

# Set source code output directory
SRC="src/chemcache"

# Set data directories
ATOMIC_INFO="${REFERENCE_DATA}/physical_data"
BASIS_SETS="${REFERENCE_DATA}/basis_sets"
DENSITIES="${REFERENCE_DATA}/atomic_densities"
MOLECULES="${REFERENCE_DATA}/molecules"

# Remove existing source files
rm -rf "${SRC}/atomic_densities"
rm -rf "${SRC}/bases"

# Recreate necessary src subdirectories
mkdir -p "${SRC}/atomic_densities"
mkdir -p "${SRC}/bases"

# Activate virtual environment
. venv/bin/activate

# Call generation scripts
echo "Calling ${REFERENCE_DATA}/generate_atomicinfo.py ${ATOMIC_INFO} ${SRC}"
${PYTHON} ${REFERENCE_DATA}/generate_atomicinfo.py ${ATOMIC_INFO} ${SRC}

# echo "Calling ${REFERENCE_DATA}/generate_densities.py"
# ${PYTHON} ${REFERENCE_DATA}/generate_densities.py ${DENSITIES} -i ${INC} ${SRC} ${TESTS} -r

echo "Calling ${REFERENCE_DATA}/generate_molecules.py ${MOLECULES} ${SRC} -r"
${PYTHON} ${REFERENCE_DATA}/generate_molecules.py ${MOLECULES} ${SRC} -r

echo "Calling ${REFERENCE_DATA}/generate_basis.py ${BASIS_SETS} ${SRC} -r"
${PYTHON} ${REFERENCE_DATA}/generate_basis.py ${BASIS_SETS} ${SRC} -r

echo "Calling ${REFERENCE_DATA}/generate_ptable_configs.py ${ATOMIC_INFO} ${SRC}"
${PYTHON} ${REFERENCE_DATA}/generate_ptable_configs.py ${ATOMIC_INFO} ${SRC}
