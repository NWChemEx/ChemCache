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

#include "cc_dash_pv_oparen_q_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_q_plus_d_cparen_z_4() {
    // Basis Set name and origin point
    std::string name("cc-pv(q+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.2000000000e-05, 7.1300000000e-04, 3.7350000000e-03,
                1.5468000000e-02, 5.2874000000e-02, 1.4569400000e-01,
                3.0268100000e-01, 4.0493600000e-01, 2.2238700000e-01,
                1.2912000000e-02, -4.6480000000e-03, 1.2280000000e-03},
      doubles_t{1.4630000000e+04, 2.1910000000e+03, 4.9820000000e+02,
                1.4090000000e+02, 4.5860000000e+01, 1.6470000000e+01,
                6.3190000000e+00, 2.5350000000e+00, 1.0350000000e+00,
                2.5280000000e-01, 1.0520000000e-01, 4.2610000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5280000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0520000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7000000000e-05, -1.3000000000e-04, -6.7900000000e-04,
                -2.8570000000e-03, -9.8130000000e-03, -2.8609000000e-02,
                -6.3760000000e-02, -1.1723100000e-01, -1.2120200000e-01,
                2.5322900000e-01, 5.8183700000e-01, 2.9195400000e-01},
      doubles_t{1.4630000000e+04, 2.1910000000e+03, 4.9820000000e+02,
                1.4090000000e+02, 4.5860000000e+01, 1.6470000000e+01,
                6.3190000000e+00, 2.5350000000e+00, 1.0350000000e+00,
                2.5280000000e-01, 1.0520000000e-01, 4.2610000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2610000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0360000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0990000000e-03, 2.5626000000e-02, 1.0376800000e-01,
                3.1154700000e-01, 4.9812800000e-01, 2.5887600000e-01},
      doubles_t{1.4030000000e+01, 3.1680000000e+00, 9.0240000000e-01,
                3.0360000000e-01, 1.1300000000e-01, 4.2860000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2860000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1288000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5690000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7680000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5210000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1310000000e-01}));
    return abs_t(name, 4, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_q_plus_d_cparen_z_4

} // namespace chemcache
