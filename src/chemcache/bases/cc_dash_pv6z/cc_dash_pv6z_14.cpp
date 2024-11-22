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

abs_t cc_dash_pv6z_14() {
    // Basis Set name and origin point
    std::string name("cc-pv6z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.7000000000e-06, 1.3600000000e-05,  7.1400000000e-05,
                2.9730000000e-04, 1.0383000000e-03,  3.1747000000e-03,
                8.7324000000e-03, 2.2383000000e-02,  5.3727300000e-02,
                1.1664900000e-01, 2.1597800000e-01,  3.0956600000e-01,
                2.8394500000e-01, 1.2223200000e-01,  1.4195200000e-02,
                3.1210000000e-04, 6.1900000000e-04,  -1.2750000000e-04,
                5.0600000000e-05, -2.7900000000e-05, 5.4000000000e-06},
      doubles_t{4.4650000000e+06, 6.6850000000e+05, 1.5220000000e+05,
                4.3300000000e+04, 1.4410000000e+04, 5.3940000000e+03,
                2.2120000000e+03, 9.6810000000e+02, 4.4120000000e+02,
                2.0710000000e+02, 9.9800000000e+01, 4.9240000000e+01,
                2.4740000000e+01, 1.2470000000e+01, 5.7950000000e+00,
                2.8300000000e+00, 1.4070000000e+00, 6.9950000000e-01,
                3.0830000000e-01, 1.3850000000e-01, 6.1450000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.0000000000e-07, -3.6000000000e-06, -1.9000000000e-05,
                -7.9100000000e-05, -2.7690000000e-04, -8.4720000000e-04,
                -2.3478000000e-03, -6.0705000000e-03, -1.4971100000e-02,
                -3.3972900000e-02, -6.9458400000e-02, -1.1900100000e-01,
                -1.5364500000e-01, -7.0468400000e-02, 2.1314900000e-01,
                4.9159600000e-01,  3.7831000000e-01,  7.5752400000e-02,
                1.5840000000e-03,  3.4350000000e-04,  5.9300000000e-05},
      doubles_t{4.4650000000e+06, 6.6850000000e+05, 1.5220000000e+05,
                4.3300000000e+04, 1.4410000000e+04, 5.3940000000e+03,
                2.2120000000e+03, 9.6810000000e+02, 4.4120000000e+02,
                2.0710000000e+02, 9.9800000000e+01, 4.9240000000e+01,
                2.4740000000e+01, 1.2470000000e+01, 5.7950000000e+00,
                2.8300000000e+00, 1.4070000000e+00, 6.9950000000e-01,
                3.0830000000e-01, 1.3850000000e-01, 6.1450000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4070000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9950000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0830000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3850000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0000000000e-07,  9.0000000000e-07,  4.9000000000e-06,
                2.0300000000e-05,  7.0900000000e-05,  2.1720000000e-04,
                6.0130000000e-04,  1.5591000000e-03,  3.8443000000e-03,
                8.7797000000e-03,  1.8038800000e-02,  3.1522400000e-02,
                4.1690500000e-02,  2.0097300000e-02,  -6.6748400000e-02,
                -1.8190600000e-01, -2.6029500000e-01, -5.4420000000e-02,
                4.2984500000e-01,  5.7876300000e-01,  1.8973500000e-01},
      doubles_t{4.4650000000e+06, 6.6850000000e+05, 1.5220000000e+05,
                4.3300000000e+04, 1.4410000000e+04, 5.3940000000e+03,
                2.2120000000e+03, 9.6810000000e+02, 4.4120000000e+02,
                2.0710000000e+02, 9.9800000000e+01, 4.9240000000e+01,
                2.4740000000e+01, 1.2470000000e+01, 5.7950000000e+00,
                2.8300000000e+00, 1.4070000000e+00, 6.9950000000e-01,
                3.0830000000e-01, 1.3850000000e-01, 6.1450000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1450000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.9900000000e-05, 5.2960000000e-04, 2.9958000000e-03,
                1.2633500000e-02, 4.1904400000e-02, 1.1025900000e-01,
                2.1883100000e-01, 3.1782800000e-01, 3.1942500000e-01,
                1.7063400000e-01, 2.6834200000e-02, -3.0200000000e-04,
                7.9320000000e-04, -1.1870000000e-04},
      doubles_t{3.5720000000e+03, 8.4600000000e+02, 2.7480000000e+02,
                1.0500000000e+02, 4.4350000000e+01, 2.0080000000e+01,
                9.5300000000e+00, 4.6340000000e+00, 2.2800000000e+00,
                1.1160000000e+00, 4.9910000000e-01, 2.2540000000e-01,
                1.0010000000e-01, 4.3320000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1160000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2540000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.2800000000e-05, -1.1260000000e-04, -6.4020000000e-04,
                -2.7029000000e-03, -9.0789000000e-03, -2.4234800000e-02,
                -4.9346000000e-02, -7.2585900000e-02, -8.0425800000e-02,
                -2.7607600000e-02, 1.7014800000e-01, 4.1115800000e-01,
                4.2507000000e-01, 1.4478900000e-01},
      doubles_t{3.5720000000e+03, 8.4600000000e+02, 2.7480000000e+02,
                1.0500000000e+02, 4.4350000000e+01, 2.0080000000e+01,
                9.5300000000e+00, 4.6340000000e+00, 2.2800000000e+00,
                1.1160000000e+00, 4.9910000000e-01, 2.2540000000e-01,
                1.0010000000e-01, 4.3320000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0010000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3320000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2386000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3767000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8530000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4880000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0580000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3510000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2250000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5750000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5280000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6310000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5150000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5570000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2310000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9460000000e-01}));
    return abs_t(name, 14, r0, shells.begin(), shells.end());
} // cc_dash_pv6z_14

} // namespace chemcache
