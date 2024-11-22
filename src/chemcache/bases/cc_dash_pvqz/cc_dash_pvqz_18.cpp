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

#include "cc_dash_pvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvqz_18() {
    // Basis Set name and origin point
    std::string name("cc-pvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.2754500000e-05, 1.7694500000e-04, 9.3128200000e-04,
                3.9286000000e-03, 1.4206400000e-02, 4.4811400000e-02,
                1.2100100000e-01, 2.6057900000e-01, 3.8136400000e-01,
                2.7605800000e-01, 5.0517900000e-02, -3.5986600000e-03,
                2.1879800000e-03, -8.6235900000e-04, 4.0522600000e-04,
                -1.1616300000e-04},
      doubles_t{9.5060000000e+05, 1.4230000000e+05, 3.2360000000e+04,
                9.1450000000e+03, 2.9700000000e+03, 1.0640000000e+03,
                4.1080000000e+02, 1.6800000000e+02, 7.1990000000e+01,
                3.1670000000e+01, 1.2890000000e+01, 5.9290000000e+00,
                2.6780000000e+00, 9.4160000000e-01, 4.2390000000e-01,
                1.7140000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.4620100000e-06, -5.0234600000e-05, -2.6480400000e-04,
                -1.1189500000e-03, -4.0827600000e-03, -1.3121600000e-02,
                -3.7285500000e-02, -8.9470900000e-02, -1.6805400000e-01,
                -1.7959400000e-01, 1.0295300000e-01, 5.6263000000e-01,
                4.5035500000e-01, 4.6060700000e-02, -6.0369100000e-03,
                1.8874400000e-03},
      doubles_t{9.5060000000e+05, 1.4230000000e+05, 3.2360000000e+04,
                9.1450000000e+03, 2.9700000000e+03, 1.0640000000e+03,
                4.1080000000e+02, 1.6800000000e+02, 7.1990000000e+01,
                3.1670000000e+01, 1.2890000000e+01, 5.9290000000e+00,
                2.6780000000e+00, 9.4160000000e-01, 4.2390000000e-01,
                1.7140000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.4160000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2390000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0205600000e-06, 1.5685100000e-05, 8.2861700000e-05,
                3.4926400000e-04, 1.2797600000e-03, 4.1036500000e-03,
                1.1778900000e-02, 2.8386800000e-02, 5.5240600000e-02,
                6.0749200000e-02, -3.6201200000e-02, -2.7539800000e-01,
                -3.6284500000e-01, 2.7311800000e-01, 6.8289900000e-01,
                2.8341700000e-01},
      doubles_t{9.5060000000e+05, 1.4230000000e+05, 3.2360000000e+04,
                9.1450000000e+03, 2.9700000000e+03, 1.0640000000e+03,
                4.1080000000e+02, 1.6800000000e+02, 7.1990000000e+01,
                3.1670000000e+01, 1.2890000000e+01, 5.9290000000e+00,
                2.6780000000e+00, 9.4160000000e-01, 4.2390000000e-01,
                1.7140000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7140000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.9575200000e-04, 4.2517200000e-03, 2.2327700000e-02,
                           8.3087800000e-02, 2.1711000000e-01, 3.7450700000e-01,
                           3.6644500000e-01, 1.2924500000e-01, 6.9224900000e-03,
                           8.3614900000e-04, -1.5266600000e-04},
                 doubles_t{1.8900000000e+03, 4.4780000000e+02, 1.4460000000e+02,
                           5.4460000000e+01, 2.2510000000e+01, 9.7740000000e+00,
                           4.3680000000e+00, 1.9590000000e+00, 8.2600000000e-01,
                           3.2970000000e-01, 1.2420000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2600000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.3886300000e-04, -1.1887000000e-03, -6.3255300000e-03,
                -2.3881300000e-02, -6.4923800000e-02, -1.1544400000e-01,
                -1.2365100000e-01, 6.4905500000e-02, 4.0363300000e-01,
                4.9059200000e-01, 2.0852800000e-01},
      doubles_t{1.8900000000e+03, 4.4780000000e+02, 1.4460000000e+02,
                5.4460000000e+01, 2.2510000000e+01, 9.7740000000e+00,
                4.3680000000e+00, 1.9590000000e+00, 8.2600000000e-01,
                3.2970000000e-01, 1.2420000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2420000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8730000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3250000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0070000000e+00}));
    return abs_t(name, 18, r0, shells.begin(), shells.end());
} // cc_dash_pvqz_18

} // namespace chemcache
