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
 * This file is autogenerated by: generate_basis_sets.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis_sets.py is run.
 */

#include "bases.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
sto-3g atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(sto_dash_3g_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(sto_dash_3g_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("sto-3g");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(1): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{3.4252509140e+00, 6.2391372980e-01, 1.6885540400e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(8): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{1.3070932140e+02, 2.3808866050e+01, 6.4436083130e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache