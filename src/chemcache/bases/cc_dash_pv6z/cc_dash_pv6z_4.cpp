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

#include "cc_dash_pv6z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv6z_4() {
    // Basis Set name and origin point
    std::string name("cc-pv6z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{6.0000000000e-06, 4.8000000000e-05, 2.5300000000e-04,
                           1.0690000000e-03, 3.8770000000e-03, 1.2459000000e-02,
                           3.5714000000e-02, 8.9575000000e-02, 1.8812500000e-01,
                           3.0929400000e-01, 3.3856900000e-01, 1.6748300000e-01,
                           1.3901000000e-02, -3.1460000000e-03,
                           1.0190000000e-03, -2.4300000000e-04},
                 doubles_t{1.2651230000e+05, 1.8942500000e+04, 4.3104100000e+03,
                           1.2208000000e+03, 3.9824000000e+02, 1.4375500000e+02,
                           5.6054600000e+01, 2.3216200000e+01, 1.0063700000e+01,
                           4.5065000000e+00, 2.0688100000e+00, 9.5398600000e-01,
                           3.4588300000e-01, 1.6745700000e-01, 7.5870000000e-02,
                           3.4438000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5398600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4588300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6745700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0000000000e-06, -9.0000000000e-06, -4.6000000000e-05,
                -1.9500000000e-04, -7.0900000000e-04, -2.2860000000e-03,
                -6.6550000000e-03, -1.7083000000e-02, -3.8114000000e-02,
                -7.0226000000e-02, -1.0876900000e-01, -1.0529500000e-01,
                9.1499000000e-02, 3.9974200000e-01, 4.8953500000e-01,
                1.5197300000e-01},
      doubles_t{1.2651230000e+05, 1.8942500000e+04, 4.3104100000e+03,
                1.2208000000e+03, 3.9824000000e+02, 1.4375500000e+02,
                5.6054600000e+01, 2.3216200000e+01, 1.0063700000e+01,
                4.5065000000e+00, 2.0688100000e+00, 9.5398600000e-01,
                3.4588300000e-01, 1.6745700000e-01, 7.5870000000e-02,
                3.4438000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.5870000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4438000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9186200000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3348000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4857700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.5800000000e-04, 2.0890000000e-03, 9.6990000000e-03,
                3.2186000000e-02, 9.4175000000e-02, 2.3165700000e-01,
                3.8527500000e-01, 3.4772300000e-01, 1.0120600000e-01},
      doubles_t{7.4085600000e+01, 1.7566700000e+01, 5.6148900000e+00,
                2.0336900000e+00, 7.9186200000e-01, 3.3348000000e-01,
                1.4857700000e-01, 6.7788000000e-02, 3.0603000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7788000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0603000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0855000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0449000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.2350000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6230000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3140000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6751000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2420000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9950000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7840000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4960000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4640000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0070000000e-01}));
    return abs_t(name, 4, r0, shells.begin(), shells.end());
} // cc_dash_pv6z_4

} // namespace chemcache
