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

abs_t mini_46() {
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
      doubles_t{6.1409097226e-02, 3.6833718336e-01, 6.8952156886e-01},
      doubles_t{1.0732268000e+04, 1.6242557000e+03, 3.5677443000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0896459293e-01, 7.2291835309e-01, 3.4452847764e-01},
      doubles_t{4.7645981000e+02, 4.6089079000e+01, 1.9011457000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.6108727603e-01, 8.4730192221e-01, 2.9803087264e-01},
      doubles_t{4.0501816000e+01, 7.0317180000e+00, 3.1102881000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.3156042504e-01, 8.0657306092e-01, 3.9062022950e-01},
      doubles_t{5.9535538000e+00, 1.2737777000e+00, 5.2932060000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.6655630422e-01, 5.9604411512e-01, 4.8897161240e-01},
      doubles_t{9.1680600000e-01, 8.0146600000e-02, 3.0334700000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{8.3539906794e-02, 4.2993413496e-01, 6.3251135144e-01},
      doubles_t{6.7616554000e+02, 1.5785238000e+02, 4.6699134000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.0359300436e-02, 4.2142830604e-01, 6.5153630935e-01},
      doubles_t{1.0125606000e+02, 1.5760330000e+01, 5.9416372000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.7277210234e-01, 5.5810970351e-01, 1.6422780103e-01},
      doubles_t{2.3082437000e+00, 9.7075110000e-01, 3.8294810000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.2512060020e-01, 4.9361430078e-01, 5.8126370092e-01},
      doubles_t{8.4456904000e+01, 2.3277977000e+01, 7.1869656000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.3265029962e-01, 5.5661549908e-01, 4.4190459927e-01},
      doubles_t{3.4886548000e+00, 1.1266213000e+00, 3.2670440000e-01}));
    return abs_t(name, 46, r0, shells.begin(), shells.end());
} // mini_46

} // namespace chemcache
