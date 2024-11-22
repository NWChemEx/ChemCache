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

abs_t aug_dash_cc_dash_pv6z_9() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv6z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{5.5600000000e-06, 4.3180000000e-05, 2.2700000000e-04,
                           9.5803000000e-04, 3.4701500000e-03, 1.1185260000e-02,
                           3.2328800000e-02, 8.2795450000e-02, 1.7988020000e-01,
                           3.0557830000e-01, 3.4026840000e-01, 1.7668240000e-01,
                           2.0854750000e-02, -7.5322000000e-04,
                           1.0744400000e-03, -7.9510000000e-05},
                 doubles_t{7.2350000000e+05, 1.0840000000e+05, 2.4680000000e+04,
                           6.9900000000e+03, 2.2820000000e+03, 8.2460000000e+02,
                           3.2180000000e+02, 1.3350000000e+02, 5.8110000000e+01,
                           2.6280000000e+01, 1.2240000000e+01, 5.7470000000e+00,
                           2.3650000000e+00, 1.0710000000e+00, 4.6810000000e-01,
                           1.9940000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7470000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3650000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0710000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.2900000000e-06, -9.9900000000e-06, -5.2600000000e-05,
                -2.2172000000e-04, -8.0692000000e-04, -2.6081700000e-03,
                -7.6740200000e-03, -2.0193530000e-02, -4.7187520000e-02,
                -9.1580090000e-02, -1.4048560000e-01, -1.0367090000e-01,
                1.5282050000e-01, 4.4614580000e-01, 4.3811810000e-01,
                1.2447050000e-01},
      doubles_t{7.2350000000e+05, 1.0840000000e+05, 2.4680000000e+04,
                6.9900000000e+03, 2.2820000000e+03, 8.2460000000e+02,
                3.2180000000e+02, 1.3350000000e+02, 5.8110000000e+01,
                2.6280000000e+01, 1.2240000000e+01, 5.7470000000e+00,
                2.3650000000e+00, 1.0710000000e+00, 4.6810000000e-01,
                1.9940000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6810000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9940000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3150000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4490000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5450000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.7721000000e-04, 1.5269100000e-03, 8.0720700000e-03,
                           3.0740210000e-02, 9.0119140000e-02, 1.9528790000e-01,
                           3.0107690000e-01, 3.3322070000e-01, 2.4114680000e-01,
                           6.9672200000e-02},
                 doubles_t{6.6000000000e+02, 1.5640000000e+02, 5.0640000000e+01,
                           1.9080000000e+01, 7.8720000000e+00, 3.4490000000e+00,
                           1.5450000000e+00, 6.8640000000e-01, 2.9860000000e-01,
                           1.2450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8640000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9860000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0573000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6130000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0130000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.5630000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3300000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4660000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7350000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7830000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1500000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.2000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0880000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9370000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5810000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7390000000e+00}));
    return abs_t(name, 9, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv6z_9

} // namespace chemcache
