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

#include "aug_dash_cc_dash_pcv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pcv5z_18() {
    // Basis Set name and origin point
    std::string name("aug-cc-pcv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{1.8000000000e-06, 1.3600000000e-05, 7.1600000000e-05,
                           3.0300000000e-04, 1.1061000000e-03, 3.6067000000e-03,
                           1.0713200000e-02, 2.9106800000e-02, 7.1660100000e-02,
                           1.5414050000e-01, 2.7041710000e-01},
                 doubles_t{7.4010000000e+06, 1.1080000000e+06, 2.5210000000e+05,
                           7.1380000000e+04, 2.3260000000e+04, 8.3900000000e+03,
                           3.2710000000e+03, 1.3570000000e+03, 5.9200000000e+02,
                           2.6910000000e+02, 1.2650000000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.0000000000e-07, -3.9000000000e-06, -2.0300000000e-05,
                -8.6100000000e-05, -3.1440000000e-04, -1.0284000000e-03,
                -3.0727000000e-03, -8.4753000000e-03, -2.1520100000e-02,
                -4.9449300000e-02, -9.8775900000e-02},
      doubles_t{7.4010000000e+06, 1.1080000000e+06, 2.5210000000e+05,
                7.1380000000e+04, 2.3260000000e+04, 8.3900000000e+03,
                3.2710000000e+03, 1.3570000000e+03, 5.9200000000e+02,
                2.6910000000e+02, 1.2650000000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1030000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9860000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4170000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0220000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5110000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7580000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3800000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.3200000000e-04, 2.0233000000e-03, 1.1034000000e-02},
      doubles_t{2.9270000000e+03, 6.9350000000e+02, 2.2470000000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5170000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5530000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5730000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1650000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3220000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4780000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7510000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6808000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8737000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9537000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8530000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7760000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9170000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5873000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5140000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.3900000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6680000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6560000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0306000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5620000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1104000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2640000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.4200000000e-01}));
    return abs_t(name, 18, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pcv5z_18

} // namespace chemcache
