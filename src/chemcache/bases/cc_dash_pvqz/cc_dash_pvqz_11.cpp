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

#include "cc_dash_pvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvqz_11() {
    // Basis Set name and origin point
    std::string name("cc-pvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.7889400000e-06, 3.7239500000e-05, 1.9583100000e-04,
                8.2669800000e-04, 3.0025100000e-03, 9.7031000000e-03,
                2.8233700000e-02, 7.3205800000e-02, 1.6289700000e-01,
                2.8870800000e-01, 3.4682900000e-01, 2.0686500000e-01,
                3.2800900000e-02, -6.4773600000e-04, 1.4587800000e-03,
                -1.7834600000e-04, 9.1478900000e-05, -8.2518200000e-05,
                2.9225200000e-05},
      doubles_t{
        1.2240000000e+06, 1.8320000000e+05, 4.1700000000e+04, 1.1810000000e+04,
        3.8530000000e+03, 1.3910000000e+03, 5.4250000000e+02, 2.2490000000e+02,
        9.7930000000e+01, 4.4310000000e+01, 2.0650000000e+01, 9.7290000000e+00,
        4.2280000000e+00, 1.9690000000e+00, 8.8900000000e-01, 3.9640000000e-01,
        6.9930000000e-02, 3.2890000000e-02, 1.6120000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1695800000e-06, -9.0911000000e-06, -4.7849900000e-05,
                -2.0196200000e-04, -7.3583700000e-04, -2.3874600000e-03,
                -7.0496900000e-03, -1.8785600000e-02, -4.4615300000e-02,
                -8.9774100000e-02, -1.4294000000e-01, -1.2431500000e-01,
                9.9964800000e-02, 4.1708000000e-01, 4.7512300000e-01,
                1.6326800000e-01, 3.0954300000e-03, -1.5526300000e-03,
                5.5776300000e-04},
      doubles_t{
        1.2240000000e+06, 1.8320000000e+05, 4.1700000000e+04, 1.1810000000e+04,
        3.8530000000e+03, 1.3910000000e+03, 5.4250000000e+02, 2.2490000000e+02,
        9.7930000000e+01, 4.4310000000e+01, 2.0650000000e+01, 9.7290000000e+00,
        4.2280000000e+00, 1.9690000000e+00, 8.8900000000e-01, 3.9640000000e-01,
        6.9930000000e-02, 3.2890000000e-02, 1.6120000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9930000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.7587100000e-07, 1.3659400000e-06, 7.1979500000e-06,
                3.0334900000e-05, 1.1075200000e-04, 3.5859600000e-04,
                1.0627200000e-03, 2.8268700000e-03, 6.7674200000e-03,
                1.3648000000e-02, 2.2281400000e-02, 1.9601100000e-02,
                -1.6770800000e-02, -7.7373400000e-02, -1.1350100000e-01,
                -1.3913000000e-01, 4.4008300000e-01, 5.3895200000e-01,
                1.3388900000e-01},
      doubles_t{
        1.2240000000e+06, 1.8320000000e+05, 4.1700000000e+04, 1.1810000000e+04,
        3.8530000000e+03, 1.3910000000e+03, 5.4250000000e+02, 2.2490000000e+02,
        9.7930000000e+01, 4.4310000000e+01, 2.0650000000e+01, 9.7290000000e+00,
        4.2280000000e+00, 1.9690000000e+00, 8.8900000000e-01, 3.9640000000e-01,
        6.9930000000e-02, 3.2890000000e-02, 1.6120000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2890000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6120000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{9.0819600000e-04, 7.4177300000e-03, 3.5746400000e-02,
                1.1852000000e-01, 2.6140300000e-01, 3.7839500000e-01,
                3.3463200000e-01, 1.2684400000e-01, -1.4711700000e-02,
                5.6865000000e-03, -1.7097400000e-03, 5.4714200000e-04},
      doubles_t{4.1340000000e+02, 9.7980000000e+01, 3.1370000000e+01,
                1.1620000000e+01, 4.6710000000e+00, 1.9180000000e+00,
                7.7750000000e-01, 3.0130000000e-01, 2.2750000000e-01,
                7.5270000000e-02, 3.1260000000e-02, 1.3420000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.5270000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1260000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-9.0174100000e-05, -7.3934200000e-04, -3.5730900000e-03,
                -1.2014200000e-02, -2.6717800000e-02, -3.9275300000e-02,
                -3.7608300000e-02, -4.3322800000e-02, 5.1800300000e-02,
                2.6019400000e-01, 5.4968100000e-01, 2.8187200000e-01},
      doubles_t{4.1340000000e+02, 9.7980000000e+01, 3.1370000000e+01,
                1.1620000000e+01, 4.6710000000e+00, 1.9180000000e+00,
                7.7750000000e-01, 3.0130000000e-01, 2.2750000000e-01,
                7.5270000000e-02, 3.1260000000e-02, 1.3420000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3420000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1160000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0390000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0060000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0090000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7760000000e-01}));
    return abs_t(name, 11, r0, shells.begin(), shells.end());
} // cc_dash_pvqz_11

} // namespace chemcache
