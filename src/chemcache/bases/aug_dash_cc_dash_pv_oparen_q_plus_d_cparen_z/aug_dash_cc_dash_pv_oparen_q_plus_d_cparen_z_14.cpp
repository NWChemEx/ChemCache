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

#include "aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_14() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(q+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.6092000000e-05, 2.0290500000e-04, 1.0671500000e-03,
                4.5059700000e-03, 1.6235900000e-02, 5.0891300000e-02,
                1.3515500000e-01, 2.8129200000e-01, 3.8533600000e-01,
                2.4565100000e-01, 3.4314500000e-02, -3.3488400000e-03,
                1.8762500000e-03, -6.9332700000e-04, 4.3830800000e-04,
                -1.2261600000e-04},
      doubles_t{5.1300000000e+05, 7.6820000000e+04, 1.7470000000e+04,
                4.9350000000e+03, 1.6020000000e+03, 5.7410000000e+02,
                2.2150000000e+02, 9.0540000000e+01, 3.8740000000e+01,
                1.6950000000e+01, 6.4520000000e+00, 2.8740000000e+00,
                1.2500000000e+00, 3.5990000000e-01, 1.6990000000e-01,
                7.0660000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.9488000000e-06, -5.3964100000e-05, -2.8471600000e-04,
                -1.2020300000e-03, -4.3839700000e-03, -1.3977600000e-02,
                -3.9351600000e-02, -9.1428300000e-02, -1.6560900000e-01,
                -1.5250500000e-01, 1.6852400000e-01, 5.6928400000e-01,
                3.9805600000e-01, 2.9150900000e-02, -8.4897600000e-03,
                1.9960300000e-03},
      doubles_t{5.1300000000e+05, 7.6820000000e+04, 1.7470000000e+04,
                4.9350000000e+03, 1.6020000000e+03, 5.7410000000e+02,
                2.2150000000e+02, 9.0540000000e+01, 3.8740000000e+01,
                1.6950000000e+01, 6.4520000000e+00, 2.8740000000e+00,
                1.2500000000e+00, 3.5990000000e-01, 1.6990000000e-01,
                7.0660000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5990000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6990000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.7806800000e-06, 1.3814800000e-05, 7.3000500000e-05,
                3.0766600000e-04, 1.1256300000e-03, 3.5843500000e-03,
                1.0172800000e-02, 2.3752000000e-02, 4.4348300000e-02,
                4.1904100000e-02, -5.0250400000e-02, -2.1657800000e-01,
                -2.8644800000e-01, 2.6325600000e-01, 6.3490000000e-01,
                2.9083200000e-01},
      doubles_t{5.1300000000e+05, 7.6820000000e+04, 1.7470000000e+04,
                4.9350000000e+03, 1.6020000000e+03, 5.7410000000e+02,
                2.2150000000e+02, 9.0540000000e+01, 3.8740000000e+01,
                1.6950000000e+01, 6.4520000000e+00, 2.8740000000e+00,
                1.2500000000e+00, 3.5990000000e-01, 1.6990000000e-01,
                7.0660000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0660000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7500000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.4814300000e-04, 3.8163900000e-03, 1.9810500000e-02,
                           7.2701700000e-02, 1.8983900000e-01, 3.3567200000e-01,
                           3.7936500000e-01, 2.0119300000e-01, 2.0851500000e-02,
                           -1.7025800000e-03, 7.5004800000e-04},
                 doubles_t{1.1220000000e+03, 2.6600000000e+02, 8.5920000000e+01,
                           3.2330000000e+01, 1.3370000000e+01, 5.8000000000e+00,
                           2.5590000000e+00, 1.1240000000e+00, 3.9880000000e-01,
                           1.5330000000e-01, 5.7280000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.9880000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5330000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-9.6488300000e-05, -8.1197100000e-04, -4.3008700000e-03,
                -1.5750200000e-02, -4.2954100000e-02, -7.5257400000e-02,
                -9.7144600000e-02, -2.2750700000e-02, 2.9198800000e-01,
                5.5067000000e-01, 2.9761800000e-01},
      doubles_t{1.1220000000e+03, 2.6600000000e+02, 8.5920000000e+01,
                3.2330000000e+01, 1.3370000000e+01, 5.8000000000e+00,
                2.5590000000e+00, 1.1240000000e+00, 3.9880000000e-01,
                1.5330000000e-01, 5.7280000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7280000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6450000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.4600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1200000000e-01}));
    return abs_t(name, 14, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_14

} // namespace chemcache
