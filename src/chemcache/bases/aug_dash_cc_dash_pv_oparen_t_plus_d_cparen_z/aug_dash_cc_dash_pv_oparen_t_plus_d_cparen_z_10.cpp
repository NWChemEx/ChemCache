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

#include "aug_dash_cc_dash_pv_oparen_t_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_t_plus_d_cparen_z_10() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(t+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{5.0200000000e-04, 3.8810000000e-03, 1.9997000000e-02,
                           7.8418000000e-02, 2.2967600000e-01, 4.3272200000e-01,
                           3.5064200000e-01, 4.3911000000e-02,
                           -7.6450000000e-03, 2.3750000000e-03},
                 doubles_t{2.4350000000e+04, 3.6500000000e+03, 8.2960000000e+02,
                           2.3400000000e+02, 7.5610000000e+01, 2.6730000000e+01,
                           9.9270000000e+00, 2.8360000000e+00, 1.1020000000e+00,
                           3.7820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8360000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1800000000e-04, -9.1500000000e-04, -4.7370000000e-03,
                -1.9233000000e-02, -6.0369000000e-02, -1.4250800000e-01,
                -1.7771000000e-01, 1.7735200000e-01, 6.0583600000e-01,
                3.6510900000e-01},
      doubles_t{2.4350000000e+04, 3.6500000000e+03, 8.2960000000e+02,
                2.3400000000e+02, 7.5610000000e+01, 2.6730000000e+01,
                9.9270000000e+00, 2.8360000000e+00, 1.1020000000e+00,
                3.7820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1330000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1430000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.7151000000e-02, 1.0765600000e-01, 3.2168100000e-01,
                           4.8523200000e-01, 3.3258400000e-01},
                 doubles_t{5.4700000000e+01, 1.2430000000e+01, 3.6790000000e+00,
                           1.1430000000e+00, 3.3000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1750000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0140000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0960000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5440000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0840000000e+00}));
    return abs_t(name, 10, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_t_plus_d_cparen_z_10

} // namespace chemcache
