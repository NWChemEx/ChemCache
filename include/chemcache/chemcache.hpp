/**
 * @file chemcache.hpp
 *
 * This is the main header of the chemcache library, defining the public API
 * of the library. This file should NOT be included in any other chemcache
 * header files, source files, or tests (except tests/chemcache.cpp, which
 * tests the functions defined in this header file).
 */

#pragma once

#include <libchemist/libchemist.hpp>

/**
 * @brief The primary namespace for the chemcache library.
 *
 * This namespace contains all functionality of the chemcache library available
 * to an end-user of the library.
 */
namespace chemcache {

/**
 * @brief Creates a basis set manager with NWX default basis sets pre-loaded
 * 
 * @return libchemist::BasisSetManager Basis set manager with defaults loaded
 */
libchemist::BasisSetManager nwx_basis_set_manager();

/**
 * @brief Creates a molecule manager with default molecules pre-loaded
 * 
 * @return libchemist::MoleculeManager Molecule manager with defaults loaded
 */
libchemist::MoleculeManager nwx_molecule_manager();

/**
 * @brief Creates a periodic table object with default elements pre-loaded
 * 
 * @return libchemist::PeriodicTable Periodic table with defaults loaded
 */
libchemist::PeriodicTable nwx_periodic_table();

} // namespace chemcache