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

abs_t aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_18() {
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
      doubles_t{2.3670000000e-04, 1.8352300000e-03, 9.5286000000e-03,
                3.8628300000e-02, 1.2408100000e-01, 2.9647100000e-01,
                4.2206800000e-01, 2.4171100000e-01, 2.0050900000e-02,
                -3.6100000000e-03, 9.7560700000e-04, -4.1131600000e-04},
      doubles_t{1.4570000000e+05, 2.1840000000e+04, 4.9720000000e+03,
                1.4080000000e+03, 4.5970000000e+02, 1.6590000000e+02,
                6.4690000000e+01, 2.6440000000e+01, 7.6280000000e+00,
                2.9960000000e+00, 6.5040000000e-01, 2.3370000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.7491000000e-05, -5.1852200000e-04, -2.7482500000e-03,
                -1.1100700000e-02, -3.8482000000e-02, -9.9759900000e-02,
                -2.0308800000e-01, -1.3560800000e-01, 5.0719500000e-01,
                6.1289800000e-01, 4.4296800000e-02, -8.9927800000e-03},
      doubles_t{1.4570000000e+05, 2.1840000000e+04, 4.9720000000e+03,
                1.4080000000e+03, 4.5970000000e+02, 1.6590000000e+02,
                6.4690000000e+01, 2.6440000000e+01, 7.6280000000e+00,
                2.9960000000e+00, 6.5040000000e-01, 2.3370000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.1045700000e-05, 1.6256500000e-04, 8.5546300000e-04,
                3.4974500000e-03, 1.2015600000e-02, 3.2136800000e-02,
                6.5527900000e-02, 4.9937000000e-02, -2.2976900000e-01,
                -4.2100600000e-01, 6.4233100000e-01, 5.6754000000e-01},
      doubles_t{1.4570000000e+05, 2.1840000000e+04, 4.9720000000e+03,
                1.4080000000e+03, 4.5970000000e+02, 1.6590000000e+02,
                6.4690000000e+01, 2.6440000000e+01, 7.6280000000e+00,
                2.9960000000e+00, 6.5040000000e-01, 2.3370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0900000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{5.7055500000e-03, 4.3046000000e-02, 1.7659100000e-01,
                           4.0686300000e-01, 4.5254900000e-01, 1.2280100000e-01,
                           -4.4599600000e-03, 2.0522500000e-03},
                 doubles_t{4.5370000000e+02, 1.0680000000e+02, 3.3730000000e+01,
                           1.2130000000e+01, 4.5940000000e+00, 1.6780000000e+00,
                           5.9090000000e-01, 1.8520000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.6065500000e-03, -1.2171400000e-02, -5.2078900000e-02,
                -1.2373700000e-01, -1.5161900000e-01, 1.4242500000e-01,
                5.8450100000e-01, 4.3754000000e-01},
      doubles_t{4.5370000000e+02, 1.0680000000e+02, 3.3730000000e+01,
                1.2130000000e+01, 4.5940000000e+00, 1.6780000000e+00,
                5.9090000000e-01, 1.8520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3300000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3900000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3800000000e-01}));
    return abs_t(name, 18, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_18

} // namespace chemcache
