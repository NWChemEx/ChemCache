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
cc-pv8z atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(cc_dash_pv8z_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(cc_dash_pv8z_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("cc-pv8z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(1): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.0805600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4762200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0757100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8968200000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{3.0100000000e-06, 5.4400000000e-06, 1.7850000000e-05,
                        7.1950000000e-05, 1.2055000000e-04, 5.2527000000e-04,
                        2.0138400000e-03, 6.8708600000e-03, 2.0954150000e-02,
                        5.6140390000e-02, 1.3038113000e-01, 2.5186989000e-01,
                        3.5225121000e-01, 2.6587036000e-01, 5.7876780000e-02},
              doubles_t{1.3940800000e+04, 3.4714300000e+03, 1.9031800000e+03,
                        6.2328600000e+02, 2.9036000000e+02, 1.2658500000e+02,
                        4.4695500000e+01, 1.6003800000e+01, 6.0805600000e+00,
                        2.4762200000e+00, 1.0757100000e+00, 4.8968200000e-01,
                        2.2978600000e-01, 1.0996550000e-01, 5.2323000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2978600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0996550000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.2323000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7369100000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9252700000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6161900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6500200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5288000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4353000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5675000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.7070200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4505600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2748800000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1628000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9436000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0380000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0609700000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.2422900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5902400000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2798500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3238000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8235900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0535100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6010500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3949000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9887600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9543500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4574300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9887600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6950300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1600000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(10): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{3.7000000000e-07, 2.7500000000e-06, 1.4090000000e-05,
                        5.8020000000e-05, 2.0791000000e-04, 6.7380000000e-04,
                        2.0046900000e-03, 5.5553600000e-03, 1.4553910000e-02,
                        3.6001190000e-02, 8.2163420000e-02, 1.6533764000e-01,
                        2.7216031000e-01, 3.1898287000e-01, 2.0602037000e-01,
                        4.6984340000e-02, 2.1556500000e-03, 1.2987500000e-03,
                        1.0780000000e-04, 8.3030000000e-05},
              doubles_t{7.7416000000e+06, 1.1997000000e+06, 2.7904900000e+05,
                        8.0647000000e+04, 2.6764900000e+04, 9.7740400000e+03,
                        3.8440400000e+03, 1.6065300000e+03, 7.0339800000e+02,
                        3.1864400000e+02, 1.4826900000e+02, 7.0616400000e+01,
                        3.4363700000e+01, 1.7055300000e+01, 8.5559700000e+00,
                        4.2083200000e+00, 2.0514000000e+00, 9.7279600000e-01,
                        4.5575300000e-01, 2.0740500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7055300000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.5559700000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2083200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0514000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.7279600000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.0000000000e-08, -6.5000000000e-07, -3.3100000000e-06,
                        -1.3620000000e-05, -4.8850000000e-05, -1.5836000000e-04,
                        -4.7203000000e-04, -1.3120100000e-03, -3.4648700000e-03,
                        -8.7041400000e-03, -2.0535580000e-02, -4.4037330000e-02,
                        -8.2457300000e-02, -1.2474347000e-01, -1.1934839000e-01,
                        2.6656460000e-02,  2.8899046000e-01,  4.4795215000e-01,
                        3.2007201000e-01,  6.6634280000e-02},
              doubles_t{7.7416000000e+06, 1.1997000000e+06, 2.7904900000e+05,
                        8.0647000000e+04, 2.6764900000e+04, 9.7740400000e+03,
                        3.8440400000e+03, 1.6065300000e+03, 7.0339800000e+02,
                        3.1864400000e+02, 1.4826900000e+02, 7.0616400000e+01,
                        3.4363700000e+01, 1.7055300000e+01, 8.5559700000e+00,
                        4.2083200000e+00, 2.0514000000e+00, 9.7279600000e-01,
                        4.5575300000e-01, 2.0740500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.5575300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0740500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5357100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7861900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9079700000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.0500000000e-06, 5.7110000000e-05, 3.3617000000e-04,
                        1.5089800000e-03, 5.5180800000e-03, 1.7059100000e-02,
                        4.5603840000e-02, 1.0301383000e-01, 1.8536807000e-01,
                        2.5916104000e-01, 2.8344532000e-01, 2.3355206000e-01,
                        1.1940731000e-01, 2.1761750000e-02},
              doubles_t{5.8653500000e+03, 1.3382600000e+03, 4.2949400000e+02,
                        1.6431400000e+02, 7.0078200000e+01, 3.2024200000e+01,
                        1.5290600000e+01, 7.5357100000e+00, 3.7861900000e+00,
                        1.9079700000e+00, 9.5277000000e-01, 4.6980400000e-01,
                        2.2854500000e-01, 1.0753300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.5277000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.6980400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2854500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0753300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4037200000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2185800000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.1777000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1318000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5876800000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.0488000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0804000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4965600000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9562600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2298500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2487500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1955200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3559000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6802500000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4815000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2812500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1610600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0908500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0877400000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5708900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8531500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4612500000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1142000000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.1819900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4100700000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.2945400000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 8,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9239400000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 9,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7880000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
