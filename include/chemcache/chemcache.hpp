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
 * @brief Load default NWChemEx basis sets into a BasisSetManager.
 *
 * @param[out] bsm BasisSetManager to load basis sets into
 */
void load_basis_sets(libchemist::BasisSetManager& bsm);

/**
 * @brief Load default NWChemEx molecules into a MoleculeManager using the
 *        given PeriodicTable instance.
 *
 * @param[out] mm MoleculeManager to load molecules into
 * @param[in]  pt PeriodicTable to retrieve atomic information from
 */
void load_molecules(libchemist::MoleculeManager& mm,
                    const libchemist::PeriodicTable& pt);

/**
 * @brief Load default NWChemEx atomic data into a PeriodicTable.
 *
 * @param[out] pt PeriodicTable to load atomic data into
 */
void load_elements(libchemist::PeriodicTable& pt);

/**
 * @brief Creates a basis set manager with NWX default basis sets pre-loaded.
 *
 * @return libchemist::BasisSetManager Basis set manager with defaults loaded
 */
libchemist::BasisSetManager nwx_basis_set_manager();

/**
 * @brief Creates a molecule manager with default molecules pre-loaded.
 *
 * @return libchemist::MoleculeManager Molecule manager with defaults loaded
 */
libchemist::MoleculeManager nwx_molecule_manager();

/**
 * @brief Creates a periodic table object with default elements pre-loaded.
 *
 * @return libchemist::PeriodicTable Periodic table with defaults loaded
 */
libchemist::PeriodicTable nwx_periodic_table();

} // namespace chemcache