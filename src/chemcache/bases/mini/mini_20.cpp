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

abs_t mini_20() {
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
      doubles_t{6.4623697912e-02, 3.7983758773e-01, 6.7832937809e-01},
      doubles_t{1.9154348000e+03, 2.8953324000e+02, 6.3106352000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.8882206403e-02, 6.6602734313e-01, 3.9991192590e-01},
      doubles_t{8.3633281000e+01, 7.5118396000e+00, 3.0145846000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.2353089302e-01, 7.0288627804e-01, 4.2355248677e-01},
      doubles_t{5.3707540000e+00, 8.3837960000e-01, 3.4622610000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7930980614e-01, 6.7287922303e-01, 4.3566641491e-01},
      doubles_t{4.2974050000e-01, 6.2173400000e-02, 2.4948000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{9.6316302765e-02, 4.4810761287e-01, 6.1992121780e-01},
      doubles_t{9.7974592000e+01, 2.2067384000e+01, 6.0938756000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.4502929208e-01, 5.8438488112e-01, 2.8766619070e-01},
      doubles_t{2.0178863000e+00, 7.6665060000e-01, 2.8431920000e-01}));
    return abs_t(name, 20, r0, shells.begin(), shells.end());
} // mini_20

} // namespace chemcache
