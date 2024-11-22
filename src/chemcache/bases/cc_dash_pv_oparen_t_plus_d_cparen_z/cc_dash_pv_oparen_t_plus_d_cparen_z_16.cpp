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

#include "cc_dash_pv_oparen_t_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_t_plus_d_cparen_z_16() {
    // Basis Set name and origin point
    std::string name("cc-pv(t+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.4214000000e-05, 4.2085500000e-04, 2.2069800000e-03,
                9.1925800000e-03, 3.2112300000e-02, 9.4668300000e-02,
                2.2363000000e-01, 3.7439300000e-01, 3.2910800000e-01,
                8.4703800000e-02, 4.4085100000e-04, 1.6482700000e-03,
                -6.2233200000e-04, 3.0130600000e-04, -8.4129000000e-05},
      doubles_t{3.7410000000e+05, 5.6050000000e+04, 1.2760000000e+04,
                3.6150000000e+03, 1.1830000000e+03, 4.2880000000e+02,
                1.6780000000e+02, 6.9470000000e+01, 2.9840000000e+01,
                1.2720000000e+01, 5.2440000000e+00, 2.2190000000e+00,
                7.7670000000e-01, 3.4900000000e-01, 1.3220000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.4983700000e-05, -1.1619800000e-04, -6.1158300000e-04,
                -2.5537000000e-03, -9.0870800000e-03, -2.7704500000e-02,
                -7.2002000000e-02, -1.4643900000e-01, -1.9515000000e-01,
                8.1919300000e-03, 5.1660100000e-01, 5.4217800000e-01,
                6.8843000000e-02, -9.1807200000e-03, 2.2683200000e-03},
      doubles_t{3.7410000000e+05, 5.6050000000e+04, 1.2760000000e+04,
                3.6150000000e+03, 1.1830000000e+03, 4.2880000000e+02,
                1.6780000000e+02, 6.9470000000e+01, 2.9840000000e+01,
                1.2720000000e+01, 5.2440000000e+00, 2.2190000000e+00,
                7.7670000000e-01, 3.4900000000e-01, 1.3220000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7670000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.3506600000e-06, 3.3714000000e-05, 1.7767400000e-04,
                7.4111600000e-04, 2.6459100000e-03, 8.0748700000e-03,
                2.1227600000e-02, 4.3832300000e-02, 6.1271600000e-02,
                -3.6151000000e-03, -2.0451000000e-01, -3.8187100000e-01,
                8.2684400000e-02, 7.1414700000e-01, 3.9379100000e-01},
      doubles_t{3.7410000000e+05, 5.6050000000e+04, 1.2760000000e+04,
                3.6150000000e+03, 1.1830000000e+03, 4.2880000000e+02,
                1.6780000000e+02, 6.9470000000e+01, 2.9840000000e+01,
                1.2720000000e+01, 5.2440000000e+00, 2.2190000000e+00,
                7.7670000000e-01, 3.4900000000e-01, 1.3220000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3220000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.4226400000e-03, 1.9279600000e-02, 8.8540100000e-02,
                2.5465400000e-01, 4.3398400000e-01, 3.5495300000e-01,
                6.1894100000e-02, -5.0297700000e-03, 2.0981300000e-03},
      doubles_t{5.7440000000e+02, 1.3580000000e+02, 4.3190000000e+01,
                1.5870000000e+01, 6.2080000000e+00, 2.4830000000e+00,
                8.6880000000e-01, 3.2290000000e-01, 1.0980000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6880000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-6.2010200000e-04, -4.9388200000e-03, -2.3264700000e-02,
                -6.8519500000e-02, -1.2389600000e-01, -9.6949900000e-02,
                2.2821500000e-01, 5.6939400000e-01, 3.6630200000e-01},
      doubles_t{5.7440000000e+02, 1.3580000000e+02, 4.3190000000e+01,
                1.5870000000e+01, 6.2080000000e+00, 2.4830000000e+00,
                8.6880000000e-01, 3.2290000000e-01, 1.0980000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0980000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7560000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.1200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5700000000e-01}));
    return abs_t(name, 16, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_t_plus_d_cparen_z_16

} // namespace chemcache
