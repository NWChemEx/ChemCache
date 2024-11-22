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

#include "cc_dash_pvdz_oparen_seg_dash_opt_cparen_.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvdz_oparen_seg_dash_opt_cparen__5() {
    // Basis Set name and origin point
    std::string name("cc-pvdz(seg-opt)");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{6.9598000000e-04, 5.3536100000e-03, 2.7134010000e-02,
                           1.0142163000e-01, 2.7219150000e-01, 4.4969078000e-01,
                           2.9205783000e-01},
                 doubles_t{4.5700000000e+03, 6.8590000000e+02, 1.5650000000e+02,
                           4.4470000000e+01, 1.4480000000e+01, 5.1310000000e+00,
                           1.8980000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.4000000000e-07, -2.3010000000e-05, -1.5757400000e-03,
                -5.1675400000e-03, -4.8767140000e-02, -7.3273470000e-02,
                5.4239928000e-01},
      doubles_t{4.5700000000e+03, 6.8590000000e+02, 4.4470000000e+01,
                1.4480000000e+01, 5.1310000000e+00, 1.8980000000e+00,
                3.3290000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0430000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.5481000000e-02, 1.9807200000e-01, 5.0523000000e-01},
      doubles_t{6.0010000000e+00, 1.2410000000e+00, 3.3640000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4300000000e-01}));
    return abs_t(name, 5, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_oparen_seg_dash_opt_cparen__5

} // namespace chemcache
