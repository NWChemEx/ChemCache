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

#include "cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvdz_3() {
    // Basis Set name and origin point
    std::string name("cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.6600000000e-04, 5.8920000000e-03, 2.9671000000e-02,
                1.0918000000e-01, 2.8278900000e-01, 4.5312300000e-01,
                2.7477400000e-01, 9.7510000000e-03, -3.1800000000e-03},
      doubles_t{1.4690000000e+03, 2.2050000000e+02, 5.0260000000e+01,
                1.4240000000e+01, 4.5810000000e+00, 1.5800000000e+00,
                5.6400000000e-01, 7.3450000000e-02, 2.8050000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.2000000000e-04, -9.2300000000e-04, -4.6890000000e-03,
                -1.7682000000e-02, -4.8902000000e-02, -9.6009000000e-02,
                -1.3638000000e-01, 5.7510200000e-01, 5.1766100000e-01},
      doubles_t{1.4690000000e+03, 2.2050000000e+02, 5.0260000000e+01,
                1.4240000000e+01, 4.5810000000e+00, 1.5800000000e+00,
                5.6400000000e-01, 7.3450000000e-02, 2.8050000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8050000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{2.2784000000e-02, 1.3910700000e-01, 5.0037500000e-01,
                           5.0847400000e-01},
                 doubles_t{1.5340000000e+00, 2.7490000000e-01, 7.3620000000e-02,
                           2.4030000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4030000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1440000000e-01}));
    return abs_t(name, 3, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_3

} // namespace chemcache
