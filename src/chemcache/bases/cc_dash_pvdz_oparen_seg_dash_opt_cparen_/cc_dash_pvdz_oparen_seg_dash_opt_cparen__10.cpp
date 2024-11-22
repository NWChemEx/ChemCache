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

abs_t cc_dash_pvdz_oparen_seg_dash_opt_cparen__10() {
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
                 doubles_t{7.3800000000e-04, 5.6778700000e-03, 2.8883000000e-02,
                           1.0860195000e-01, 2.9112758000e-01, 4.5018937000e-01,
                           2.5926297000e-01},
                 doubles_t{1.7880000000e+04, 2.6830000000e+03, 6.1150000000e+02,
                           1.7350000000e+02, 5.6640000000e+01, 2.0420000000e+01,
                           7.8100000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.4000000000e-07, -3.2830000000e-05, -2.3458900000e-03,
                -8.3535400000e-03, -7.0654800000e-02, -4.6853080000e-02,
                5.7056347000e-01},
      doubles_t{1.7880000000e+04, 2.6830000000e+03, 1.7350000000e+02,
                5.6640000000e+01, 2.0420000000e+01, 7.8100000000e+00,
                1.6530000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8690000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.6087000000e-02, 2.4018100000e-01, 5.0874400000e-01},
      doubles_t{2.8390000000e+01, 6.2700000000e+00, 1.6950000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3170000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2020000000e+00}));
    return abs_t(name, 10, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_oparen_seg_dash_opt_cparen__10

} // namespace chemcache
