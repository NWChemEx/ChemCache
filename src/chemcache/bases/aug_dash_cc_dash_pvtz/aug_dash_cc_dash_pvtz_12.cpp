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

abs_t aug_dash_cc_dash_pvtz_12() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.2992900000e-05, 5.6665200000e-04, 2.9626900000e-03,
                1.2296200000e-02, 4.2732400000e-02, 1.2301300000e-01,
                2.7483200000e-01, 4.0181800000e-01, 2.6469700000e-01,
                3.3261200000e-02, -4.4133500000e-03, 2.0602400000e-03,
                -8.3847000000e-04, 7.0819500000e-04, -2.0060200000e-04},
      doubles_t{1.6490000000e+05, 2.4710000000e+04, 5.6280000000e+03,
                1.5960000000e+03, 5.2100000000e+02, 1.8800000000e+02,
                7.3010000000e+01, 2.9900000000e+01, 1.2540000000e+01,
                4.3060000000e+00, 1.8260000000e+00, 7.4170000000e-01,
                1.4570000000e-01, 7.6120000000e-02, 3.3100000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.8424800000e-05, -1.4350000000e-04, -7.4871000000e-04,
                -3.1440700000e-03, -1.1048100000e-02, -3.3605800000e-02,
                -8.2594600000e-02, -1.5931400000e-01, -1.5288800000e-01,
                1.9084900000e-01, 5.7996400000e-01, 3.7202900000e-01,
                2.1169100000e-02, -1.1934400000e-02, 2.7732900000e-03},
      doubles_t{1.6490000000e+05, 2.4710000000e+04, 5.6280000000e+03,
                1.5960000000e+03, 5.2100000000e+02, 1.8800000000e+02,
                7.3010000000e+01, 2.9900000000e+01, 1.2540000000e+01,
                4.3060000000e+00, 1.8260000000e+00, 7.4170000000e-01,
                1.4570000000e-01, 7.6120000000e-02, 3.3100000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4570000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{3.5517600000e-06, 2.7642000000e-05, 1.4440400000e-04,
                6.0574400000e-04, 2.1352700000e-03, 6.4993400000e-03,
                1.6144600000e-02, 3.1576600000e-02, 3.1637400000e-02,
                -4.3914000000e-02, -1.5109300000e-01, -2.1766800000e-01,
                2.5980700000e-01, 5.4724500000e-01, 3.2591000000e-01},
      doubles_t{1.6490000000e+05, 2.4710000000e+04, 5.6280000000e+03,
                1.5960000000e+03, 5.2100000000e+02, 1.8800000000e+02,
                7.3010000000e+01, 2.9900000000e+01, 1.2540000000e+01,
                4.3060000000e+00, 1.8260000000e+00, 7.4170000000e-01,
                1.4570000000e-01, 7.6120000000e-02, 3.3100000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3100000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2000000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{2.0753200000e-03, 1.6286900000e-02, 7.3869700000e-02,
                           2.1429700000e-01, 3.8215400000e-01, 3.9817800000e-01,
                           1.5287800000e-01, -4.3754000000e-03,
                           6.2958200000e-03, -1.5676900000e-03},
                 doubles_t{3.1690000000e+02, 7.4860000000e+01, 2.3720000000e+01,
                           8.6690000000e+00, 3.3630000000e+00, 1.3100000000e+00,
                           4.9110000000e-01, 2.3640000000e-01, 8.7330000000e-02,
                           3.2370000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7330000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.2972700000e-04, -2.5875400000e-03, -1.1912000000e-02,
                -3.5022700000e-02, -6.3996800000e-02, -7.0443600000e-02,
                -3.7583600000e-02, 1.7704300000e-01, 5.6670700000e-01,
                3.9614300000e-01},
      doubles_t{3.1690000000e+02, 7.4860000000e+01, 2.3720000000e+01,
                8.6690000000e+00, 3.3630000000e+00, 1.3100000000e+00,
                4.9110000000e-01, 2.3640000000e-01, 8.7330000000e-02,
                3.2370000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2370000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.0000000000e-03}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9980000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4170000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6500000000e-02}));
    return abs_t(name, 12, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvtz_12

} // namespace chemcache
