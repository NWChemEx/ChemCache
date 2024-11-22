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

#include "cc_dash_pwcvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pwcvtz_6() {
    // Basis Set name and origin point
    std::string name("cc-pwcvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{8.8280000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{5.3100000000e-04, 4.1080000000e-03, 2.1087000000e-02,
                           8.1853000000e-02, 2.3481700000e-01, 4.3440100000e-01,
                           3.4612900000e-01, 3.9378000000e-02,
                           -8.9830000000e-03, 2.3850000000e-03},
                 doubles_t{8.2360000000e+03, 1.2350000000e+03, 2.8080000000e+02,
                           7.9270000000e+01, 2.5590000000e+01, 8.9970000000e+00,
                           3.3190000000e+00, 9.0590000000e-01, 3.6430000000e-01,
                           1.2850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.0590000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1300000000e-04, -8.7800000000e-04, -4.5400000000e-03,
                -1.8133000000e-02, -5.5760000000e-02, -1.2689500000e-01,
                -1.7035200000e-01, 1.4038200000e-01, 5.9868400000e-01,
                3.9538900000e-01},
      doubles_t{8.2360000000e+03, 1.2350000000e+03, 2.8080000000e+02,
                7.9270000000e+01, 2.5590000000e+01, 8.9970000000e+00,
                3.3190000000e+00, 9.0590000000e-01, 3.6430000000e-01,
                1.2850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4740000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1552000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4420000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8270000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.4031000000e-02, 8.6866000000e-02, 2.9021600000e-01,
                           5.0100800000e-01, 3.4340600000e-01},
                 doubles_t{1.8710000000e+01, 4.1330000000e+00, 1.2000000000e+00,
                           3.8270000000e-01, 1.2090000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2090000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9500000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0970000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6100000000e-01}));
    return abs_t(name, 6, r0, shells.begin(), shells.end());
} // cc_dash_pwcvtz_6

} // namespace chemcache
