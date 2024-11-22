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

#include "aug_dash_cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvdz_9() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.2100000000e-04, 5.5530000000e-03, 2.8267000000e-02,
                1.0644400000e-01, 2.8681400000e-01, 4.4864100000e-01,
                2.6476100000e-01, 1.5333000000e-02, -2.3320000000e-03},
      doubles_t{1.4710000000e+04, 2.2070000000e+03, 5.0280000000e+02,
                1.4260000000e+02, 4.6470000000e+01, 1.6700000000e+01,
                6.3560000000e+00, 1.3160000000e+00, 3.8970000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.6500000000e-04, -1.3080000000e-03, -6.4950000000e-03,
                -2.6691000000e-02, -7.3690000000e-02, -1.7077600000e-01,
                -1.1232700000e-01, 5.6281400000e-01, 5.6877800000e-01},
      doubles_t{1.4710000000e+04, 2.2070000000e+03, 5.0280000000e+02,
                1.4260000000e+02, 4.6470000000e+01, 1.6700000000e+01,
                6.3560000000e+00, 1.3160000000e+00, 3.8970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8630000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.4878000000e-02, 2.3571800000e-01, 5.0852100000e-01,
                           4.5812000000e-01},
                 doubles_t{2.2670000000e+01, 4.9770000000e+00, 1.3470000000e+00,
                           3.4710000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4710000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5020000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6400000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6400000000e-01}));
    return abs_t(name, 9, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvdz_9

} // namespace chemcache
