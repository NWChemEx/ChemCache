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

#include "aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_8() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(d+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.1000000000e-04, 5.4700000000e-03, 2.7837000000e-02,
                1.0480000000e-01, 2.8306200000e-01, 4.4871900000e-01,
                2.7095200000e-01, 1.5458000000e-02, -2.5850000000e-03},
      doubles_t{1.1720000000e+04, 1.7590000000e+03, 4.0080000000e+02,
                1.1370000000e+02, 3.7030000000e+01, 1.3270000000e+01,
                5.0250000000e+00, 1.0130000000e+00, 3.0230000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.6000000000e-04, -1.2630000000e-03, -6.2670000000e-03,
                -2.5716000000e-02, -7.0924000000e-02, -1.6541100000e-01,
                -1.1695500000e-01, 5.5736800000e-01, 5.7275900000e-01},
      doubles_t{1.1720000000e+04, 1.7590000000e+03, 4.0080000000e+02,
                1.1370000000e+02, 3.7030000000e+01, 1.3270000000e+01,
                5.0250000000e+00, 1.0130000000e+00, 3.0230000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0230000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8960000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.3018000000e-02, 2.2891300000e-01, 5.0872800000e-01,
                           4.6053100000e-01},
                 doubles_t{1.7700000000e+01, 3.8540000000e+00, 1.0460000000e+00,
                           2.7530000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7530000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8560000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1850000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3200000000e-01}));
    return abs_t(name, 8, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_8

} // namespace chemcache
