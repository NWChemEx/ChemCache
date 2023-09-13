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

#include "bases.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using molecular_basis_pt = simde::MolecularBasisSet;
using molecular_basis_t  = simde::type::ao_basis_set;
using atomic_basis_pt    = simde::AtomicBasisSetFromZ;
using atomic_basis_t     = simde::type::atomic_basis_set;

static constexpr auto module_desc = R"(
Molecular Basis Set
---------------------------------

This module returns a molecular basis set for the input molecule. The basis set
is determined by the "Atomic Basis" submodule.
)";

MODULE_CTOR(molecular_basis) {
    description(module_desc);
    satisfies_property_type<molecular_basis_pt>();
    add_submodule<atomic_basis_pt>("Atomic Basis");
}

MODULE_RUN(molecular_basis) {
    const auto& [mol] = molecular_basis_pt::unwrap_inputs(inputs);
    auto& atoms_mod   = submods.at("Atomic Basis");

    atomic_basis_t ci;
    molecular_basis_t aobs;
    for(const auto& ai : mol) {
        ci = atoms_mod.run_as<atomic_basis_pt>(ai.Z());
        for(auto i : {0, 1, 2}) ci.center().coord(i) = ai.coord(i);
        aobs.add_center(ci);
    }

    auto rv = results();
    return molecular_basis_pt::wrap_results(rv, aobs);
}

} // namespace chemcache
