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

abs_t cc_dash_pvdz_6() {
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
      doubles_t{6.9200000000e-04, 5.3290000000e-03, 2.7077000000e-02,
                1.0171800000e-01, 2.7474000000e-01, 4.4856400000e-01,
                2.8507400000e-01, 1.5204000000e-02, -3.1910000000e-03},
      doubles_t{6.6650000000e+03, 1.0000000000e+03, 2.2800000000e+02,
                6.4710000000e+01, 2.1060000000e+01, 7.4950000000e+00,
                2.7970000000e+00, 5.2150000000e-01, 1.5960000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.4600000000e-04, -1.1540000000e-03, -5.7250000000e-03,
                -2.3312000000e-02, -6.3955000000e-02, -1.4998100000e-01,
                -1.2726200000e-01, 5.4452900000e-01, 5.8049600000e-01},
      doubles_t{6.6650000000e+03, 1.0000000000e+03, 2.2800000000e+02,
                6.4710000000e+01, 2.1060000000e+01, 7.4950000000e+00,
                2.7970000000e+00, 5.2150000000e-01, 1.5960000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5960000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{3.8109000000e-02, 2.0948000000e-01, 5.0855700000e-01,
                           4.6884200000e-01},
                 doubles_t{9.4390000000e+00, 2.0020000000e+00, 5.4560000000e-01,
                           1.5170000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5170000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5000000000e-01}));
    return abs_t(name, 6, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_6

} // namespace chemcache
