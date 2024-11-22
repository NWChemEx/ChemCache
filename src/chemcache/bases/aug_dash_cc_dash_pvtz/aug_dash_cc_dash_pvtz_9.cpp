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

#include "aug_dash_cc_dash_pvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvtz_9() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{5.0700000000e-04, 3.9230000000e-03, 2.0200000000e-02,
                           7.9010000000e-02, 2.3043900000e-01, 4.3287200000e-01,
                           3.4996400000e-01, 4.3233000000e-02,
                           -7.8920000000e-03, 2.3840000000e-03},
                 doubles_t{1.9500000000e+04, 2.9230000000e+03, 6.6450000000e+02,
                           1.8750000000e+02, 6.0620000000e+01, 2.1420000000e+01,
                           7.9500000000e+00, 2.2570000000e+00, 8.8150000000e-01,
                           3.0410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2570000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1700000000e-04, -9.1200000000e-04, -4.7170000000e-03,
                -1.9086000000e-02, -5.9655000000e-02, -1.4001000000e-01,
                -1.7678200000e-01, 1.7162500000e-01, 6.0504300000e-01,
                3.6951200000e-01},
      doubles_t{1.9500000000e+04, 2.9230000000e+03, 6.6450000000e+02,
                1.8750000000e+02, 6.0620000000e+01, 2.1420000000e+01,
                7.9500000000e+00, 2.2570000000e+00, 8.8150000000e-01,
                3.0410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1580000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1320000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.6665000000e-02, 1.0447200000e-01, 3.1726000000e-01,
                           4.8734300000e-01, 3.3460400000e-01},
                 doubles_t{4.3880000000e+01, 9.9260000000e+00, 2.9300000000e+00,
                           9.1320000000e-01, 2.6720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3610000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1070000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9170000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2400000000e-01}));
    return abs_t(name, 9, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvtz_9

} // namespace chemcache
