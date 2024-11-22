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

#include "cc_dash_pvtz_oparen_seg_dash_opt_cparen_.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvtz_oparen_seg_dash_opt_cparen__15() {
    // Basis Set name and origin point
    std::string name("cc-pvtz(seg-opt)");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{3.4145540000e-02, 1.0020529000e-01, 2.3440516000e-01,
                           3.8261335000e-01, 3.1852212000e-01, 7.0529890000e-02,
                           -4.1259000000e-03},
                 doubles_t{9.8680000000e+02, 3.5740000000e+02, 1.3960000000e+02,
                           5.7630000000e+01, 2.4600000000e+01, 1.0120000000e+01,
                           4.2830000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.6737000000e-04, -1.5776500000e-03, -1.0679790000e-02,
                -4.6427290000e-02, -1.0482865000e-01, 5.6928250000e-02,
                5.3986949000e-01, 5.2891840000e-01},
      doubles_t{9.8680000000e+02, 3.5740000000e+02, 1.3960000000e+02,
                5.7630000000e+01, 2.4600000000e+01, 1.0120000000e+01,
                4.2830000000e+00, 1.8050000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1580000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.1860000000e-05, -3.4880000000e-04, -1.7392100000e-03,
                -2.1049700000e-03, 6.1789000000e-04, -2.5815620000e-02,
                -2.2812942000e-01, 8.8113942000e-01},
      doubles_t{3.5740000000e+02, 1.3960000000e+02, 5.7630000000e+01,
                2.4600000000e+01, 1.0120000000e+01, 4.2830000000e+00,
                1.8050000000e+00, 2.7820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0550000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.3849800000e-03, 1.8919050000e-02, 8.6700490000e-02,
                2.4958297000e-01, 4.3133755000e-01, 3.7583743000e-01},
      doubles_t{5.0490000000e+02, 1.1940000000e+02, 3.7960000000e+01,
                1.3950000000e+01, 5.4570000000e+00, 2.1770000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0100000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-7.1570000000e-05, -5.8994000000e-04, -4.7976400000e-03,
                -1.1134600000e-02, -2.3694960000e-02, 7.5524298000e-01},
      doubles_t{1.1940000000e+02, 3.7960000000e+01, 1.3950000000e+01,
                5.4570000000e+00, 2.1770000000e+00, 2.8770000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.7140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5200000000e-01}));
    return abs_t(name, 15, r0, shells.begin(), shells.end());
} // cc_dash_pvtz_oparen_seg_dash_opt_cparen__15

} // namespace chemcache
