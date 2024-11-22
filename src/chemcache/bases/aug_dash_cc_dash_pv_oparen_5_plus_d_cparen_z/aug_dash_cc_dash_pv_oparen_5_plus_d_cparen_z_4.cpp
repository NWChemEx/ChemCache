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

#include "aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z_4() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(5+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{1.8000000000e-05, 1.3800000000e-04, 7.2300000000e-04,
                           3.0390000000e-03, 1.0908000000e-02, 3.4035000000e-02,
                           9.1193000000e-02, 1.9926800000e-01, 3.2935500000e-01,
                           3.4048900000e-01, 1.4374000000e-01, 6.5170000000e-03,
                           -1.8650000000e-03, 4.8000000000e-04},
                 doubles_t{5.4620000000e+04, 8.1800000000e+03, 1.8620000000e+03,
                           5.2730000000e+02, 1.7200000000e+02, 6.2100000000e+01,
                           2.4210000000e+01, 9.9930000000e+00, 4.3050000000e+00,
                           1.9210000000e+00, 8.6630000000e-01, 2.4750000000e-01,
                           1.0090000000e-01, 4.1290000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6630000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4750000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0090000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.0000000000e-06, -2.5000000000e-05, -1.3100000000e-04,
                -5.5800000000e-04, -1.9880000000e-03, -6.3700000000e-03,
                -1.7217000000e-02, -4.0858000000e-02, -7.4237000000e-02,
                -1.1923400000e-01, -8.7777000000e-02, 2.7878800000e-01,
                5.8349900000e-01, 2.6900700000e-01},
      doubles_t{5.4620000000e+04, 8.1800000000e+03, 1.8620000000e+03,
                5.2730000000e+02, 1.7200000000e+02, 6.2100000000e+01,
                2.4210000000e+01, 9.9930000000e+00, 4.3050000000e+00,
                1.9210000000e+00, 8.6630000000e-01, 2.4750000000e-01,
                1.0090000000e-01, 4.1290000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1290000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3340000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8080000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{6.3300000000e-04, 4.8080000000e-03, 2.0527000000e-02,
                           6.7816000000e-02, 1.9286800000e-01, 3.7869200000e-01,
                           4.0057400000e-01, 1.3806700000e-01},
                 doubles_t{4.3750000000e+01, 1.0330000000e+01, 3.2260000000e+00,
                           1.1270000000e+00, 4.3340000000e-01, 1.8080000000e-01,
                           7.8270000000e-02, 3.3720000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8270000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3720000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6000000000e-03}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7175000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6460000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4040000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5150000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1270000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5680000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0780000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1830000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{9.7900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1420000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0330000000e-01}));
    return abs_t(name, 4, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z_4

} // namespace chemcache
