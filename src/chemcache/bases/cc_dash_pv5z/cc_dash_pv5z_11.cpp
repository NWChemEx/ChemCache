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

#include "cc_dash_pv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv5z_11() {
    // Basis Set name and origin point
    std::string name("cc-pv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0000000000e-06,  1.8000000000e-05,  9.5000000000e-05,
                4.0100000000e-04,  1.4590000000e-03,  4.7460000000e-03,
                1.4031000000e-02,  3.7733000000e-02,  9.0702000000e-02,
                1.8646600000e-01,  3.0183700000e-01,  3.2383000000e-01,
                1.6870000000e-01,  2.3083000000e-02,  -1.4700000000e-04,
                1.0990000000e-03,  -1.0100000000e-04, 6.1000000000e-05,
                -5.3000000000e-05, 1.8000000000e-05},
      doubles_t{2.1855720000e+06, 3.2722840000e+05, 7.4466840000e+04,
                2.1093150000e+04, 6.8818980000e+03, 2.4846960000e+03,
                9.6922320000e+02, 4.0206430000e+02, 1.7535450000e+02,
                7.9651990000e+01, 3.7386720000e+01, 1.8001940000e+01,
                8.7243710000e+00, 3.8577150000e+00, 1.8156860000e+00,
                8.3825400000e-01, 3.8193500000e-01, 7.1679000000e-02,
                3.3916000000e-02, 1.6525000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0000000000e-06, -4.0000000000e-06, -2.3000000000e-05,
                -9.8000000000e-05, -3.5700000000e-04, -1.1650000000e-03,
                -3.4640000000e-03, -9.4950000000e-03, -2.3587000000e-02,
                -5.2394000000e-02, -9.8028000000e-02, -1.4367300000e-01,
                -1.0229800000e-01, 1.3802900000e-01,  4.2900600000e-01,
                4.4781300000e-01,  1.4454300000e-01,  2.6170000000e-03,
                -1.1920000000e-03, 4.4000000000e-04},
      doubles_t{2.1855720000e+06, 3.2722840000e+05, 7.4466840000e+04,
                2.1093150000e+04, 6.8818980000e+03, 2.4846960000e+03,
                9.6922320000e+02, 4.0206430000e+02, 1.7535450000e+02,
                7.9651990000e+01, 3.7386720000e+01, 1.8001940000e+01,
                8.7243710000e+00, 3.8577150000e+00, 1.8156860000e+00,
                8.3825400000e-01, 3.8193500000e-01, 7.1679000000e-02,
                3.3916000000e-02, 1.6525000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8193500000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1679000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0000000000e-06, 3.0000000000e-06, 1.5000000000e-05,
                5.4000000000e-05, 1.7500000000e-04, 5.2000000000e-04,
                1.4310000000e-03, 3.5540000000e-03, 7.9540000000e-03,
                1.4955000000e-02, 2.2472000000e-02, 1.6205000000e-02,
                -2.3550000000e-02, -8.0120000000e-02, -1.1292100000e-01,
                -1.3298900000e-01, 4.1751600000e-01, 5.4944300000e-01,
                1.4811000000e-01},
      doubles_t{
        3.2722840000e+05, 7.4466840000e+04, 2.1093150000e+04, 6.8818980000e+03,
        2.4846960000e+03, 9.6922320000e+02, 4.0206430000e+02, 1.7535450000e+02,
        7.9651990000e+01, 3.7386720000e+01, 1.8001940000e+01, 8.7243710000e+00,
        3.8577150000e+00, 1.8156860000e+00, 8.3825400000e-01, 3.8193500000e-01,
        7.1679000000e-02, 3.3916000000e-02, 1.6525000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3916000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6525000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.6200000000e-04, 1.4080000000e-03, 7.5860000000e-03,
                2.9615000000e-02, 8.8477000000e-02, 1.9551900000e-01,
                3.0662100000e-01, 3.4124200000e-01, 2.2969000000e-01,
                5.4921000000e-02, -3.4410000000e-03, 2.3400000000e-03,
                -6.6200000000e-04, 1.9900000000e-04},
      doubles_t{1.1195780000e+03, 2.6532390000e+02, 8.5995530000e+01,
                3.2537590000e+01, 1.3515650000e+01, 5.9668560000e+00,
                2.7000450000e+00, 1.2185120000e+00, 5.4218700000e-01,
                2.2741300000e-01, 1.3304000000e-01, 5.7577000000e-02,
                2.5971000000e-02, 1.1901000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3304000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7577000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.6000000000e-05, -1.4000000000e-04, -7.5400000000e-04,
                -2.9680000000e-03, -8.9190000000e-03, -2.0029000000e-02,
                -3.1380000000e-02, -3.6489000000e-02, -2.9625000000e-02,
                -1.0146000000e-02, 8.2128000000e-02, 3.3455300000e-01,
                5.0814200000e-01, 1.9663600000e-01},
      doubles_t{1.1195780000e+03, 2.6532390000e+02, 8.5995530000e+01,
                3.2537590000e+01, 1.3515650000e+01, 5.9668560000e+00,
                2.7000450000e+00, 1.2185120000e+00, 5.4218700000e-01,
                2.2741300000e-01, 1.3304000000e-01, 5.7577000000e-02,
                2.5971000000e-02, 1.1901000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5971000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1901000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8980000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2830000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0640000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4930000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3570000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3180000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2220000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2400000000e-01}));
    return abs_t(name, 11, r0, shells.begin(), shells.end());
} // cc_dash_pv5z_11

} // namespace chemcache
