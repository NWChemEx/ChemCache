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

#include "d_dash_aug_dash_cc_dash_pv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t d_dash_aug_dash_cc_dash_pv5z_10() {
    // Basis Set name and origin point
    std::string name("d-aug-cc-pv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{2.6000000000e-05, 2.0000000000e-04, 1.0500000000e-03,
                           4.4000000000e-03, 1.5649000000e-02, 4.7758000000e-02,
                           1.2294300000e-01, 2.5248300000e-01, 3.6631400000e-01,
                           2.7961700000e-01, 6.1651000000e-02, 9.3400000000e-04,
                           1.3670000000e-03, -3.8000000000e-05},
                 doubles_t{2.6270000000e+05, 3.9350000000e+04, 8.9550000000e+03,
                           2.5380000000e+03, 8.2990000000e+02, 3.0150000000e+02,
                           1.1900000000e+02, 5.0000000000e+01, 2.1980000000e+01,
                           9.8910000000e+00, 4.3270000000e+00, 1.8040000000e+00,
                           7.2880000000e-01, 2.8670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3270000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8040000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.0000000000e-06, -4.7000000000e-05, -2.4700000000e-04,
                -1.0380000000e-03, -3.7110000000e-03, -1.1593000000e-02,
                -3.1086000000e-02, -7.0972000000e-02, -1.2726600000e-01,
                -1.5123100000e-01, 2.4666000000e-02, 3.9961200000e-01,
                5.2661600000e-01, 1.9085400000e-01},
      doubles_t{2.6270000000e+05, 3.9350000000e+04, 8.9550000000e+03,
                2.5380000000e+03, 8.2990000000e+02, 3.0150000000e+02,
                1.1900000000e+02, 5.0000000000e+01, 2.1980000000e+01,
                9.8910000000e+00, 4.3270000000e+00, 1.8040000000e+00,
                7.2880000000e-01, 2.8670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2880000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2690000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3150000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.0380000000e-03, 8.3750000000e-03, 3.9693000000e-02,
                           1.2805600000e-01, 2.7032500000e-01, 3.7096300000e-01,
                           3.2849900000e-01, 1.2758800000e-01},
                 doubles_t{2.9910000000e+02, 7.0730000000e+01, 2.2480000000e+01,
                           8.2460000000e+00, 3.2690000000e+00, 1.3150000000e+00,
                           5.1580000000e-01, 1.9180000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1580000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9180000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2300000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8370000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8440000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5020000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7300000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0900000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7380000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0570000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4600000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8800000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7760000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6280000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0200000000e-01}));
    return abs_t(name, 10, r0, shells.begin(), shells.end());
} // d_dash_aug_dash_cc_dash_pv5z_10

} // namespace chemcache
