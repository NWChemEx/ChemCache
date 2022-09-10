/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file chemcache.hpp
 *
 * This is the main header of the chemcache library, defining the public API
 * of the library. This file should NOT be included in any other chemcache
 * header files, source files, or tests (except tests/chemcache.cpp, which
 * tests the functions defined in this header file).
 */

#pragma once

#include <chemist/chemist.hpp>

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
void load_basis_sets(chemist::BasisSetManager& bsm);

/**
 * @brief Load default NWChemEx molecules into a MoleculeManager using the
 *        given PeriodicTable instance.
 *
 * @param[out] mm MoleculeManager to load molecules into
 * @param[in]  pt PeriodicTable to retrieve atomic information from
 */
void load_molecules(chemist::MoleculeManager& mm,
                    const chemist::PeriodicTable& pt);

/**
 * @brief Load default NWChemEx atomic data into a PeriodicTable.
 *
 * @param[out] pt PeriodicTable to load atomic data into
 */
void load_elements(chemist::PeriodicTable& pt);

/**
 * @brief Creates a basis set manager with NWX default basis sets pre-loaded.
 *
 * @return chemist::BasisSetManager Basis set manager with defaults loaded
 */
chemist::BasisSetManager nwx_basis_set_manager();

/**
 * @brief Creates a molecule manager with default molecules pre-loaded.
 *
 * @return chemist::MoleculeManager Molecule manager with defaults loaded
 */
chemist::MoleculeManager nwx_molecule_manager();

/**
 * @brief Creates a periodic table object with default elements pre-loaded.
 *
 * @return chemist::PeriodicTable Periodic table with defaults loaded
 */
chemist::PeriodicTable nwx_periodic_table();

} // namespace chemcache