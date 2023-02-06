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
 * @file chemcache.cpp
 *
 * @brief Implementation of chemcache.hpp API
 */

#include "chemcache/chemcache.hpp"

namespace chemcache {

chemist::BasisSetManager nwx_basis_set_manager() {
    chemist::BasisSetManager bsm;
    load_basis_sets(bsm);

    return bsm;
}

chemist::MoleculeManager nwx_molecule_manager() {
    chemist::MoleculeManager mm;
    chemist::PeriodicTable pt;

    load_elements(pt);
    load_molecules(mm, pt);

    return mm;
}

chemist::PeriodicTable nwx_periodic_table() {
    chemist::PeriodicTable pt;
    load_elements(pt);
    load_elec_configs(pt);

    return pt;
}

} // namespace chemcache
