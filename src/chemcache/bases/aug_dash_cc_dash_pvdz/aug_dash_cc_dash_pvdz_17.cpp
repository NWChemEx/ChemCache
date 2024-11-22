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

#include "aug_dash_cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvdz_17() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4115300000e-04, 1.8709500000e-03, 9.7082700000e-03,
                3.9315300000e-02, 1.2593200000e-01, 2.9934100000e-01,
                4.2188600000e-01, 2.3720100000e-01, 1.9153100000e-02,
                -3.3479200000e-03, 9.2988300000e-04, -3.9637900000e-04},
      doubles_t{1.2790000000e+05, 1.9170000000e+04, 4.3630000000e+03,
                1.2360000000e+03, 4.0360000000e+02, 1.4570000000e+02,
                5.6810000000e+01, 2.3230000000e+01, 6.6440000000e+00,
                2.5750000000e+00, 5.3710000000e-01, 1.9380000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.7892200000e-05, -5.2183600000e-04, -2.7651300000e-03,
                -1.1153700000e-02, -3.8591900000e-02, -9.9484800000e-02,
                -2.0139200000e-01, -1.3031300000e-01, 5.0944300000e-01,
                6.1072500000e-01, 4.2154900000e-02, -9.2342700000e-03},
      doubles_t{1.2790000000e+05, 1.9170000000e+04, 4.3630000000e+03,
                1.2360000000e+03, 4.0360000000e+02, 1.4570000000e+02,
                5.6810000000e+01, 2.3230000000e+01, 6.6440000000e+00,
                2.5750000000e+00, 5.3710000000e-01, 1.9380000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0498600000e-05, 1.5829800000e-04, 8.3363900000e-04,
                3.3988000000e-03, 1.1673800000e-02, 3.0962200000e-02,
                6.2953300000e-02, 4.6025700000e-02, -2.1931200000e-01,
                -4.0877300000e-01, 6.3846500000e-01, 5.6236200000e-01},
      doubles_t{1.2790000000e+05, 1.9170000000e+04, 4.3630000000e+03,
                1.2360000000e+03, 4.0360000000e+02, 1.4570000000e+02,
                5.6810000000e+01, 2.3230000000e+01, 6.6440000000e+00,
                2.5750000000e+00, 5.3710000000e-01, 1.9380000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9380000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0800000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{5.2598200000e-03, 3.9833200000e-02, 1.6465500000e-01,
                           3.8732200000e-01, 4.5707200000e-01, 1.5163600000e-01,
                           1.8161500000e-03, 1.8829600000e-03},
                 doubles_t{4.1760000000e+02, 9.8330000000e+01, 3.1040000000e+01,
                           1.1190000000e+01, 4.2490000000e+00, 1.6240000000e+00,
                           5.3220000000e-01, 1.6200000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.4357000000e-03, -1.0779600000e-02, -4.7007500000e-02,
                -1.1103000000e-01, -1.5327500000e-01, 8.9460900000e-02,
                5.7944400000e-01, 4.8327200000e-01},
      doubles_t{4.1760000000e+02, 9.8330000000e+01, 3.1040000000e+01,
                1.1190000000e+01, 4.2490000000e+00, 1.6240000000e+00,
                5.3220000000e-01, 1.6200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9600000000e-01}));
    return abs_t(name, 17, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvdz_17

} // namespace chemcache
