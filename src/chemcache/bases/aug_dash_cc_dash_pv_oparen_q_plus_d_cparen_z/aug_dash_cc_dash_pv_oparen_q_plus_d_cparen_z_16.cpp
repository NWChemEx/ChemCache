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

abs_t aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_16() {
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
      doubles_t{2.3602500000e-05, 1.8348200000e-04, 9.6427800000e-04,
                4.0653700000e-03, 1.4697300000e-02, 4.6508100000e-02,
                1.2550800000e-01, 2.6843300000e-01, 3.8480900000e-01,
                2.6537200000e-01, 4.3732600000e-02, -3.7880700000e-03,
                2.1808300000e-03, -8.3694400000e-04, 4.4808800000e-04,
                -1.2529300000e-04},
      doubles_t{7.2780000000e+05, 1.0900000000e+05, 2.4800000000e+04,
                7.0140000000e+03, 2.2780000000e+03, 8.1470000000e+02,
                3.1340000000e+02, 1.2770000000e+02, 5.4480000000e+01,
                2.3850000000e+01, 9.4280000000e+00, 4.2900000000e+00,
                1.9090000000e+00, 6.2700000000e-01, 2.8730000000e-01,
                1.1720000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.5217900000e-06, -5.0663100000e-05, -2.6683300000e-04,
                -1.1260100000e-03, -4.1118600000e-03, -1.3245400000e-02,
                -3.7700400000e-02, -8.9855400000e-02, -1.6709800000e-01,
                -1.6935400000e-01, 1.2782400000e-01, 5.6486200000e-01,
                4.3176700000e-01, 3.8939800000e-02, -7.3026000000e-03,
                1.9232700000e-03},
      doubles_t{7.2780000000e+05, 1.0900000000e+05, 2.4800000000e+04,
                7.0140000000e+03, 2.2780000000e+03, 8.1470000000e+02,
                3.1340000000e+02, 1.2770000000e+02, 5.4480000000e+01,
                2.3850000000e+01, 9.4280000000e+00, 4.2900000000e+00,
                1.9090000000e+00, 6.2700000000e-01, 2.8730000000e-01,
                1.1720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8730000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.8940600000e-06, 1.4694800000e-05, 7.7546000000e-05,
                3.2650900000e-04, 1.1968600000e-03, 3.8479900000e-03,
                1.1053900000e-02, 2.6464500000e-02, 5.0877100000e-02,
                5.3003000000e-02, -4.2551800000e-02, -2.5085300000e-01,
                -3.3315200000e-01, 2.6379600000e-01, 6.6684900000e-01,
                2.8845100000e-01},
      doubles_t{7.2780000000e+05, 1.0900000000e+05, 2.4800000000e+04,
                7.0140000000e+03, 2.2780000000e+03, 8.1470000000e+02,
                3.1340000000e+02, 1.2770000000e+02, 5.4480000000e+01,
                2.3850000000e+01, 9.4280000000e+00, 4.2900000000e+00,
                1.9090000000e+00, 6.2700000000e-01, 2.8730000000e-01,
                1.1720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1720000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2800000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.4118300000e-04, 3.7757100000e-03, 1.9836000000e-02,
                           7.4206300000e-02, 1.9732700000e-01, 3.5185100000e-01,
                           3.7868700000e-01, 1.7093100000e-01, 1.5158700000e-02,
                           6.7198100000e-05, 4.0549100000e-04},
                 doubles_t{1.5460000000e+03, 3.6640000000e+02, 1.1840000000e+02,
                           4.4530000000e+01, 1.8380000000e+01, 7.9650000000e+00,
                           3.5410000000e+00, 1.5910000000e+00, 6.2050000000e-01,
                           2.4200000000e-01, 9.0140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2050000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.1311000000e-04, -9.5858100000e-04, -5.1347100000e-03,
                -1.9264100000e-02, -5.3598000000e-02, -9.6033300000e-02,
                -1.1818300000e-01, 9.2319400000e-03, 3.5884100000e-01,
                5.2581800000e-01, 2.4887200000e-01},
      doubles_t{1.5460000000e+03, 3.6640000000e+02, 1.1840000000e+02,
                4.4530000000e+01, 1.8380000000e+01, 7.9650000000e+00,
                3.5410000000e+00, 1.5910000000e+00, 6.2050000000e-01,
                2.4200000000e-01, 9.0140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.0140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1590000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0190000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9700000000e-01}));
    return abs_t(name, 16, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_16

} // namespace chemcache
