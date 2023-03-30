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
REF_DATA="reference_data"

# Set source code output directories
SRC="src/chemcache"
ATOM_DEN="${SRC}/atomic_densities"
ATOM="${SRC}/atoms"
BASES="${SRC}/bases"
ELEC_CONFIGS="${SRC}/electronic_configurations"
MOLES="${SRC}/molecules"

# Set data directories
ATOMIC_INFO="${REF_DATA}/physical_data"
BASIS_SETS="${REF_DATA}/basis_sets"
DENSITIES="${REF_DATA}/atomic_densities/default"
MOLECULES="${REF_DATA}/molecules"
TESTS="tests/chemcache"

# Remove existing source files
rm -rf $ATOM_DEN
rm -rf $ATOM
rm -rf $BASES
rm -rf $ELEC_CONFIGS
rm -rf $MOLES

# Recreate necessary src subdirectories
mkdir -p $ATOM_DEN
mkdir -p $ATOM
mkdir -p $BASES
mkdir -p $ELEC_CONFIGS
mkdir -p $MOLES

# Activate virtual environment
. venv/bin/activate

# Call generation scripts
echo "Calling ${REF_DATA}/generate_atomicinfo.py ${ATOMIC_INFO} ${ATOM}"
${PYTHON} ${REF_DATA}/generate_atomicinfo.py ${ATOMIC_INFO} ${ATOM}

echo "Calling ${REF_DATA}/generate_densities.py ${DENSITIES} ${ATOM_DEN} ${TESTS} -r" 
${PYTHON} ${REF_DATA}/generate_densities.py ${DENSITIES} ${ATOM_DEN} ${TESTS} -r

echo "Calling ${REF_DATA}/generate_molecules.py ${MOLECULES} ${MOLES} -r"
${PYTHON} ${REF_DATA}/generate_molecules.py ${MOLECULES} ${MOLES} -r

echo "Calling ${REF_DATA}/generate_basis.py ${BASIS_SETS} ${BASES} -r"
${PYTHON} ${REF_DATA}/generate_basis.py ${BASIS_SETS} ${BASES} -r

echo "Calling ${REF_DATA}/generate_elec_configs.py ${ATOMIC_INFO} ${ELEC_CONFIGS}"
${PYTHON} ${REF_DATA}/generate_elec_configs.py ${ATOMIC_INFO} ${ELEC_CONFIGS}
