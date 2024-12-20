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

#include "sto_dash_3g.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t sto_dash_3g_12() {
    // Basis Set name and origin point
    std::string name("sto-3g");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
      doubles_t{2.9923741370e+02, 5.4506468450e+01, 1.4751577520e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
      doubles_t{1.5121823520e+01, 3.5139865790e+00, 1.1428574980e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
      doubles_t{1.5121823520e+01, 3.5139865790e+00, 1.1428574980e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
      doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
      doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
    return abs_t(name, 12, r0, shells.begin(), shells.end());
} // sto_dash_3g_12

} // namespace chemcache
