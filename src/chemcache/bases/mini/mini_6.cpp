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

abs_t mini_6() {
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
      doubles_t{7.0740101715e-02, 3.9538040958e-01, 6.6331071608e-01},
      doubles_t{1.5317226000e+02, 2.3073030000e+01, 4.9232900000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-8.1380500552e-02, 5.7485320390e-01, 5.0241280341e-01},
      doubles_t{5.7255700000e+00, 4.5504000000e-01, 1.4707000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0993060168e-01, 4.6271270705e-01, 6.2751370956e-01},
      doubles_t{4.2513100000e+00, 8.6327000000e-01, 2.0135000000e-01}));
    return abs_t(name, 6, r0, shells.begin(), shells.end());
} // mini_6

} // namespace chemcache
