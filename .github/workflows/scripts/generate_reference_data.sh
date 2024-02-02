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

# Set data directories
REF_DATA="reference_data"
ATOMIC_INFO="${REF_DATA}/physical_data"
BASIS_SETS="${REF_DATA}/basis_sets"
DENSITIES="${REF_DATA}/atomic_densities/default"
MOLECULES="${REF_DATA}/molecules"

# Set source code output directories
SRC="src/chemcache"
ATOM="${SRC}/atoms"
BASES="${SRC}/bases"
MOLES="${SRC}/molecules"
EXP_SRC="experimental/src/chemcache"
ATOM_DEN="${EXP_SRC}/atomic_densities"
ELEC_CONFIGS="${EXP_SRC}/electronic_configurations"

# Activate virtual environment
. venv/bin/activate

# Call generation script(s)
SCRIPTS_DIR="utils/data_management"

echo "Calling ${SCRIPTS_DIR}/generate_atomicinfo.py ${ATOMIC_INFO} ${ATOM}"
python ${SCRIPTS_DIR}/generate_atomicinfo.py ${ATOMIC_INFO} ${ATOM}

echo "Calling ${SCRIPTS_DIR}/generate_basis.py ${BASIS_SETS} ${BASES} -r"
python ${SCRIPTS_DIR}/generate_basis.py ${BASIS_SETS} ${BASES} -r

echo "Calling ${SCRIPTS_DIR}/generate_molecules.py ${MOLECULES} ${MOLES} -r"
python ${SCRIPTS_DIR}/generate_molecules.py ${MOLECULES} ${MOLES} -r

echo "Calling ${SCRIPTS_DIR}/generate_densities.py ${DENSITIES} ${ATOM_DEN} -r" 
python ${SCRIPTS_DIR}/generate_densities.py ${DENSITIES} ${ATOM_DEN} -r

echo "Calling ${SCRIPTS_DIR}/generate_elec_configs.py ${ATOMIC_INFO} ${ELEC_CONFIGS}"
python ${SCRIPTS_DIR}/generate_elec_configs.py ${ATOMIC_INFO} ${ELEC_CONFIGS}
