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
 * This file is autogenerated by: generate_basis.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis.py is run.
 */

#include "bases.hpp"
#include <simde/basis_sets/atomic_basis_set.hpp>
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
5-21g atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(five_dash_21g_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(five_dash_21g_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("5-21g");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(1): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0, doubles_t{1.5628497870e-01, 9.0469087670e-01},
              doubles_t{5.4471780000e+00, 8.2454724000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8319158000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(2): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0, doubles_t{1.7522987180e-01, 8.9348234650e-01},
              doubles_t{1.3626700000e+01, 1.9993500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.8299300000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(3): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{6.1218469100e-03, 4.5112961500e-02, 1.9269415000e-01,
                        4.6854420800e-01, 4.4060751500e-01},
              doubles_t{2.7539444400e+02, 4.1435175400e+01, 9.3669937800e+00,
                        2.5377253300e+00, 7.4663654000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0, doubles_t{-2.5253679580e-01, 1.0973407950e+00},
              doubles_t{7.3456426060e-01, 8.7197964050e-02}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1, doubles_t{1.4359173210e-01, 9.4780305060e-01},
              doubles_t{7.3456426060e-01, 8.7197964050e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0438657410e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0438657410e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(4): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{5.4099743720e-03, 4.0251532530e-02, 1.7685814290e-01,
                        4.5255936570e-01, 4.7029338000e-01},
              doubles_t{5.5401000000e+02, 8.3263100000e+01, 1.8863500000e+01,
                        5.1778200000e+00, 1.5560200000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0, doubles_t{-4.7742905160e-01, 1.2474501350e+00},
              doubles_t{1.4417524910e+00, 3.0186105970e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1, doubles_t{2.0114189920e-01, 8.8448255690e-01},
              doubles_t{1.4417524910e+00, 3.0186105970e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0096138750e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0096138750e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
