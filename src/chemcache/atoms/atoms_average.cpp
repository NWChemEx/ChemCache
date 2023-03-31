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

/*
 * This file is autogenerated by: generate_atomicinfo.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_atomicinfo.py is run.
 */

#include "atoms.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using atom_pt = simde::AtomFromZ;
using atom_t  = simde::type::atom;

static constexpr auto module_desc = R"(
Atoms with Abundance-Weighted Mass
---------------------------------

This module returns an atom with abundance-weighted mass.
This module was autogenerated.
)";

MODULE_CTOR(atoms_average) {
    description(module_desc);
    satisfies_property_type<atom_pt>();
}

MODULE_RUN(atoms_average) {
    const auto& [Z] = atom_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {
        case(1): {
            atom_t atom{"H", 1ul, 1837.4260218693814, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(8): {
            atom_t atom{"O", 8ul, 29165.122045980286, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        default: {
            throw std::out_of_range("Atom not available for Z");
        }
    }
}

} // namespace chemcache
