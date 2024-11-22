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

#include "aug_dash_cc_dash_pv6z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv6z_13() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv6z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9000000000e-06, 1.4500000000e-05,  7.6200000000e-05,
                3.1580000000e-04, 1.0974000000e-03,  3.3697000000e-03,
                9.3222000000e-03, 2.3799200000e-02,  5.6819100000e-02,
                1.2246800000e-01, 2.2389700000e-01,  3.1344600000e-01,
                2.7497500000e-01, 1.1056400000e-01,  1.1921500000e-02,
                6.5280000000e-04, 4.7350000000e-04,  -2.7000000000e-05,
                1.7600000000e-05, -1.0300000000e-05, 1.9000000000e-06},
      doubles_t{3.6520000000e+06, 5.4680000000e+05, 1.2450000000e+05,
                3.5440000000e+04, 1.1840000000e+04, 4.4340000000e+03,
                1.8120000000e+03, 7.9150000000e+02, 3.6100000000e+02,
                1.6950000000e+02, 8.1680000000e+01, 4.0280000000e+01,
                2.0250000000e+01, 1.0230000000e+01, 4.8020000000e+00,
                2.3390000000e+00, 1.1630000000e+00, 5.8820000000e-01,
                2.3110000000e-01, 1.0270000000e-01, 4.5210000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.0000000000e-07, -3.8000000000e-06, -1.9800000000e-05,
                -8.2100000000e-05, -2.8580000000e-04, -8.7850000000e-04,
                -2.4482000000e-03, -6.3100000000e-03, -1.5485400000e-02,
                -3.4958900000e-02, -7.0772900000e-02, -1.1942300000e-01,
                -1.4884200000e-01, -5.9046500000e-02, 2.1669300000e-01,
                4.7655700000e-01,  3.7559000000e-01,  8.6658100000e-02,
                2.4337000000e-03,  -2.3910000000e-04, 1.1980000000e-04},
      doubles_t{3.6520000000e+06, 5.4680000000e+05, 1.2450000000e+05,
                3.5440000000e+04, 1.1840000000e+04, 4.4340000000e+03,
                1.8120000000e+03, 7.9150000000e+02, 3.6100000000e+02,
                1.6950000000e+02, 8.1680000000e+01, 4.0280000000e+01,
                2.0250000000e+01, 1.0230000000e+01, 4.8020000000e+00,
                2.3390000000e+00, 1.1630000000e+00, 5.8820000000e-01,
                2.3110000000e-01, 1.0270000000e-01, 4.5210000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1630000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3110000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0270000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0000000000e-07,  9.0000000000e-07,  4.6000000000e-06,
                1.9000000000e-05,  6.5900000000e-05,  2.0310000000e-04,
                5.6470000000e-04,  1.4620000000e-03,  3.5794000000e-03,
                8.1516000000e-03,  1.6527600000e-02,  2.8546700000e-02,
                3.6148400000e-02,  1.5380400000e-02,  -6.1214100000e-02,
                -1.5126300000e-01, -2.2301200000e-01, -8.3286400000e-02,
                3.8683000000e-01,  5.9043200000e-01,  2.1007500000e-01},
      doubles_t{3.6520000000e+06, 5.4680000000e+05, 1.2450000000e+05,
                3.5440000000e+04, 1.1840000000e+04, 4.4340000000e+03,
                1.8120000000e+03, 7.9150000000e+02, 3.6100000000e+02,
                1.6950000000e+02, 8.1680000000e+01, 4.0280000000e+01,
                2.0250000000e+01, 1.0230000000e+01, 4.8020000000e+00,
                2.3390000000e+00, 1.1630000000e+00, 5.8820000000e-01,
                2.3110000000e-01, 1.0270000000e-01, 4.5210000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5210000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7370000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{6.3800000000e-05, 5.6310000000e-04, 3.1691000000e-03,
                           1.3240100000e-02, 4.3340300000e-02, 1.1195000000e-01,
                           2.1779600000e-01, 3.1167500000e-01, 3.1672200000e-01,
                           1.7961300000e-01, 3.1510900000e-02, 7.2800000000e-05,
                           6.5930000000e-04, -1.4140000000e-04},
                 doubles_t{2.8840000000e+03, 6.8320000000e+02, 2.2200000000e+02,
                           8.4820000000e+01, 3.5810000000e+01, 1.6220000000e+01,
                           7.7020000000e+00, 3.7410000000e+00, 1.8310000000e+00,
                           8.8780000000e-01, 3.9890000000e-01, 1.7180000000e-01,
                           7.2980000000e-02, 3.0690000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.8780000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.9890000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7180000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2980000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-8.0000000000e-06, -6.5100000000e-05, -3.9990000000e-04,
                -1.5369000000e-03, -5.5644000000e-03, -1.3110600000e-02,
                -2.9720000000e-02, -3.4719500000e-02, -5.5162100000e-02,
                -2.1699000000e-03, 2.5156200000e-02, 3.1722700000e-01,
                9.7780300000e-02, 7.1634400000e-01},
      doubles_t{2.8840000000e+03, 6.8320000000e+02, 2.2200000000e+02,
                8.4820000000e+01, 3.5810000000e+01, 1.6220000000e+01,
                7.7020000000e+00, 3.7410000000e+00, 1.8310000000e+00,
                8.8780000000e-01, 3.9890000000e-01, 1.7180000000e-01,
                7.2980000000e-02, 3.0690000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0690000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0210000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2143000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.4490000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0320000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7210000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3430000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6660000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7560000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2840000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6250000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7710000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0460000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5450000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3020000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9900000000e-01}));
    return abs_t(name, 13, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv6z_13

} // namespace chemcache
