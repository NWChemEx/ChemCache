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

abs_t cc_dash_pv_oparen_q_plus_d_cparen_z_3() {
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
      doubles_t{1.1700000000e-04, 9.1100000000e-04, 4.7280000000e-03,
                1.9197000000e-02, 6.3047000000e-02, 1.6320800000e-01,
                3.1482700000e-01, 3.9393600000e-01, 1.9691800000e-01,
                9.9970000000e-03, -5.4020000000e-03, 1.7040000000e-03},
      doubles_t{6.6010000000e+03, 9.8970000000e+02, 2.2570000000e+02,
                6.4290000000e+01, 2.1180000000e+01, 7.7240000000e+00,
                3.0030000000e+00, 1.2120000000e+00, 4.9300000000e-01,
                9.5150000000e-02, 4.7910000000e-02, 2.2200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5150000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7910000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.8000000000e-05, -1.4200000000e-04, -7.4100000000e-04,
                -3.0200000000e-03, -1.0123000000e-02, -2.7094000000e-02,
                -5.7359000000e-02, -9.3895000000e-02, -1.2109100000e-01,
                2.7660800000e-01, 5.4954800000e-01, 2.7738500000e-01},
      doubles_t{6.6010000000e+03, 9.8970000000e+02, 2.2570000000e+02,
                6.4290000000e+01, 2.1180000000e+01, 7.7240000000e+00,
                3.0030000000e+00, 1.2120000000e+00, 4.9300000000e-01,
                9.5150000000e-02, 4.7910000000e-02, 2.2200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1920000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.3880000000e-03, 1.9316000000e-02, 7.9104000000e-02,
                2.7409500000e-01, 5.1948800000e-01, 2.8442300000e-01},
      doubles_t{6.2500000000e+00, 1.3700000000e+00, 3.6720000000e-01,
                1.1920000000e-01, 4.4740000000e-02, 1.7950000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4740000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7950000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1350000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3190000000e-01}));
    return abs_t(name, 3, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_q_plus_d_cparen_z_3

} // namespace chemcache
