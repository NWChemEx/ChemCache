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

abs_t mini_66() {
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
      doubles_t{6.0349699564e-02, 3.6420079737e-01, 6.9355589499e-01},
      doubles_t{2.2512038000e+04, 3.4087516000e+03, 7.5021737000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1164461036e-01, 7.2846096759e-01, 3.3823783138e-01},
      doubles_t{9.9666912000e+02, 9.8727277000e+01, 4.2083165000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9223789568e-01, 8.5650658733e-01, 3.1130169540e-01},
      doubles_t{8.9147986000e+01, 1.7526082000e+01, 8.1030013000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.8836417473e-01, 7.0412755418e-01, 5.3320906530e-01},
      doubles_t{1.5711507000e+01, 4.2363528000e+00, 2.2096177000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0063900115e-02, 2.1642990247e-01, 8.0161350915e-01},
      doubles_t{4.8826712000e+00, 7.7030950000e-01, 3.9223240000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.8449200633e-02, 2.6555230214e-01, 7.4302750600e-01},
      doubles_t{9.2097320000e-01, 7.1245100000e-02, 2.7101100000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{7.7244396692e-02, 4.1510538222e-01, 6.4741017227e-01},
      doubles_t{1.5272025000e+03, 3.5881603000e+02, 1.0760264000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.5458598290e-02, 3.8948328121e-01, 6.7691156735e-01},
      doubles_t{2.0687090000e+02, 4.0194535000e+01, 1.6221668000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.2346488935e-01, 5.3576298236e-01, 2.1360969297e-01},
      doubles_t{7.8825557000e+00, 3.7978777000e+00, 1.8177367000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.5718999646e-02, 5.4489338772e-01, 5.1923498830e-01},
      doubles_t{2.1190102000e+00, 8.0074650000e-01, 2.8941950000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.0650780339e-01, 4.7360601506e-01, 5.9588561894e-01},
      doubles_t{2.3167736000e+02, 6.6102671000e+01, 2.1916854000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.8453730477e-01, 5.7442700963e-01, 3.0198030506e-01},
      doubles_t{1.2545087000e+01, 4.7409523000e+00, 1.8318939000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3,
                 doubles_t{8.2445596452e-02, 3.2497248602e-01, 5.2360247747e-01,
                           4.2687078163e-01},
                 doubles_t{4.2262676000e+01, 1.2869755000e+01, 4.2445315000e+00,
                           1.2421218000e+00}));
    return abs_t(name, 66, r0, shells.begin(), shells.end());
} // mini_66

} // namespace chemcache
