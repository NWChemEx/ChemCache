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

abs_t sto_dash_3g_37(){

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
              doubles_t{2.9386015290e+03, 5.3526993680e+02, 1.4486493420e+02}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{2.4850703690e+02, 5.7747691050e+01, 1.8781341420e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{2.4850703690e+02, 5.7747691050e+01, 1.8781341420e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.2776350230e-01, 2.1754360440e-01, 9.1667696110e-01},
              doubles_t{2.3505340970e+01, 7.1698782010e+00, 2.7663619090e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.9515111550e-03, 5.7776646910e-01, 4.8464603660e-01},
              doubles_t{2.3505340970e+01, 7.1698782010e+00, 2.7663619090e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-3.0884412140e-01, 1.9606411650e-02, 1.1310344420e+00},
              doubles_t{2.2477968200e+00, 8.2957839350e-01, 3.6635056530e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.2154686000e-01, 5.7152276040e-01, 5.4989494710e-01},
              doubles_t{2.2477968200e+00, 8.2957839350e-01, 3.6635056530e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-3.8426426080e-01, -1.9725674380e-01, 1.3754955120e+00},
              doubles_t{4.8699399190e-01, 2.6221615650e-01, 1.1582548750e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-3.4816915260e-01, 6.2903236900e-01, 6.6628327430e-01},
              doubles_t{4.8699399190e-01, 2.6221615650e-01, 1.1582548750e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 2,
              doubles_t{2.1976795080e-01, 6.5554736270e-01, 2.8657325900e-01},
              doubles_t{2.3505340970e+01, 7.1698782010e+00, 2.7663619090e+00}));
     return abs_t(name, 37, r0, shells.begin(), shells.end());
} // sto_dash_3g_37

} // chemcache
