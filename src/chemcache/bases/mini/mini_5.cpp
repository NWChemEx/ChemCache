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

#include "mini.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t mini_5() {
    // Basis Set name and origin point
    std::string name("mini");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.8650500316e-02, 3.8993290179e-01, 6.7139510309e-01},
      doubles_t{1.0843704000e+02, 1.6120560000e+01, 3.3734300000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-8.2419195645e-02, 5.5906437046e-01, 5.1679477269e-01},
      doubles_t{3.5665000000e+00, 2.9547000000e-01, 9.8050000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0590010835e-01, 4.5718053605e-01, 6.3186074982e-01},
      doubles_t{2.5720700000e+00, 5.1694000000e-01, 1.2314000000e-01}));
    return abs_t(name, 5, r0, shells.begin(), shells.end());
} // mini_5

} // namespace chemcache
