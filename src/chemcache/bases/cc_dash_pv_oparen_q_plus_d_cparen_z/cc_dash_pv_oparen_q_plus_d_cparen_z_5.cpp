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

abs_t cc_dash_pv_oparen_q_plus_d_cparen_z_5() {
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
      doubles_t{8.8000000000e-05, 6.8700000000e-04, 3.6000000000e-03,
                1.4949000000e-02, 5.1435000000e-02, 1.4330200000e-01,
                3.0093500000e-01, 4.0352600000e-01, 2.2534000000e-01,
                1.5407000000e-02, -3.9550000000e-03, 1.1240000000e-03},
      doubles_t{2.3870000000e+04, 3.5750000000e+03, 8.1280000000e+02,
                2.2970000000e+02, 7.4690000000e+01, 2.6810000000e+01,
                1.0320000000e+01, 4.1780000000e+00, 1.7270000000e+00,
                4.7040000000e-01, 1.8960000000e-01, 7.3940000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7040000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8960000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.8000000000e-05, -1.3900000000e-04, -7.2500000000e-04,
                -3.0630000000e-03, -1.0581000000e-02, -3.1365000000e-02,
                -7.1012000000e-02, -1.3210300000e-01, -1.2307200000e-01,
                2.6181900000e-01, 5.8666200000e-01, 2.9049400000e-01},
      doubles_t{2.3870000000e+04, 3.5750000000e+03, 8.1280000000e+02,
                2.2970000000e+02, 7.4690000000e+01, 2.6810000000e+01,
                1.0320000000e+01, 4.1780000000e+00, 1.7270000000e+00,
                4.7040000000e-01, 1.8960000000e-01, 7.3940000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3940000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0710000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.0950000000e-03, 3.3206000000e-02, 1.3231400000e-01,
                3.3181800000e-01, 4.7206300000e-01, 2.5797900000e-01},
      doubles_t{2.2260000000e+01, 5.0580000000e+00, 1.4870000000e+00,
                5.0710000000e-01, 1.8120000000e-01, 6.4630000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8120000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4630000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1100000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.8200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7300000000e-01}));
    return abs_t(name, 5, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_q_plus_d_cparen_z_5

} // namespace chemcache
