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
#include <simde/simde.hpp>

namespace chemcache {

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
cc-pv9z atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(cc_dash_pv9z_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(cc_dash_pv9z_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("cc-pv9z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(10): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{7.6022733000e-07, 4.1194508000e-06, 1.8538108000e-05,
                        7.1219188000e-05, 2.3862683000e-04, 7.2101492000e-04,
                        2.0419452000e-03, 5.5616242000e-03, 1.4639531000e-02,
                        3.6543265000e-02, 8.3450105000e-02, 1.6635182300e-01,
                        2.7004955000e-01, 3.1356768000e-01, 2.0546336700e-01,
                        5.1063387000e-02, 2.9358374000e-03, 1.3729280000e-03,
                        1.2730329000e-04, 1.5282880000e-04, -1.1012535000e-05},
              doubles_t{3.4064780000e+06, 7.6640500000e+05, 2.0948000000e+05,
                        6.5557000000e+04, 2.2993000000e+04, 8.8550800000e+03,
                        3.6509900000e+03, 1.5732200000e+03, 6.9700500000e+02,
                        3.1565200000e+02, 1.4644300000e+02, 6.9787400000e+01,
                        3.4150100000e+01, 1.7119100000e+01, 8.7124390000e+00,
                        4.3897130000e+00, 2.1972220000e+00, 1.0830370000e+00,
                        5.3398700000e-01, 2.6287000000e-01, 1.1903900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7119100000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.7124390000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3897130000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1972220000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0830370000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.7853772000e-07, -9.6732804000e-07, -4.3539656000e-06,
                        -1.6725280000e-05, -5.6068386000e-05, -1.6947290000e-04,
                        -4.8085599000e-04, -1.3136901000e-03, -3.4855166000e-03,
                        -8.8381618000e-03, -2.0871370000e-02, -4.4373979000e-02,
                        -8.1942970000e-02, -1.2239692100e-01, -1.1801170000e-01,
                        1.3627065000e-02,  2.5378522600e-01,  4.1560987100e-01,
                        3.4253565100e-01,  1.1538413400e-01,  6.0878761000e-03},
              doubles_t{3.4064780000e+06, 7.6640500000e+05, 2.0948000000e+05,
                        6.5557000000e+04, 2.2993000000e+04, 8.8550800000e+03,
                        3.6509900000e+03, 1.5732200000e+03, 6.9700500000e+02,
                        3.1565200000e+02, 1.4644300000e+02, 6.9787400000e+01,
                        3.4150100000e+01, 1.7119100000e+01, 8.7124390000e+00,
                        4.3897130000e+00, 2.1972220000e+00, 1.0830370000e+00,
                        5.3398700000e-01, 2.6287000000e-01, 1.1903900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3398700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6287000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1903900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5090010000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7621380000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8983760000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.7507069000e-06, 1.3632768000e-05, 7.5339633000e-05,
                        3.6253157000e-04, 1.5198830000e-03, 5.4886444000e-03,
                        1.7020535000e-02, 4.5765537000e-02, 1.0390945600e-01,
                        1.8659152200e-01, 2.5763062700e-01, 2.7731142900e-01,
                        2.2826900400e-01, 1.2289081400e-01, 2.8353858000e-02,
                        6.2584132000e-04},
              doubles_t{1.1543600000e+04, 2.9121500000e+03, 9.9768400000e+02,
                        3.8364300000e+02, 1.5853100000e+02, 6.9534200000e+01,
                        3.2006000000e+01, 1.5280500000e+01, 7.5090010000e+00,
                        3.7621380000e+00, 1.8983760000e+00, 9.5757980000e-01,
                        4.8207580000e-01, 2.4197130000e-01, 1.1957590000e-01,
                        4.8533100000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.5757980000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8207580000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4197130000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1957590000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8533100000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4589270000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8055740000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.4251750000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9199810000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5682510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3406380000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.9981900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6530800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9567000000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0781160000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9393720000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2720160000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8025620000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.9303600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.4706600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1368360000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1436460000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.1208540000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2759130000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7532860000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.3836800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4430230000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.8682310000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2902340000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3392940000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2755240000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4820150000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6470650000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9458160000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0360050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0584000000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.0400000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4000000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 9,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.4532250000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 9,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4606900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 10,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4900000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
