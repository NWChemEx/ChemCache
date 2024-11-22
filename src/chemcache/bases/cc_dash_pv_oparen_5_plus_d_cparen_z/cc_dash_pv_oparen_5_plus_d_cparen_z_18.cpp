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

#include "cc_dash_pv_oparen_5_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_5_plus_d_cparen_z_18() {
    // Basis Set name and origin point
    std::string name("cc-pv(5+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.7501800000e-06,  1.3614700000e-05,  7.1629700000e-05,
                3.0303500000e-04,  1.1061000000e-03,  3.6067900000e-03,
                1.0713400000e-02,  2.9107400000e-02,  7.1661700000e-02,
                1.5414400000e-01,  2.7042300000e-01,  3.3486200000e-01,
                2.2434700000e-01,  5.0008100000e-02,  1.4972600000e-04,
                2.1036900000e-03,  -4.1072600000e-04, -8.8861700000e-04,
                -1.5623600000e-03, -4.2683800000e-04},
      doubles_t{7.4010000000e+06, 1.1080000000e+06, 2.5210000000e+05,
                7.1380000000e+04, 2.3260000000e+04, 8.3900000000e+03,
                3.2710000000e+03, 1.3570000000e+03, 5.9200000000e+02,
                2.6910000000e+02, 1.2650000000e+02, 6.1030000000e+01,
                2.9860000000e+01, 1.4170000000e+01, 7.0220000000e+00,
                3.5110000000e+00, 1.7580000000e+00, 7.8410000000e-01,
                3.4800000000e-01, 1.4910000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.2044800000e-07, -4.0490900000e-06, -2.1302300000e-05,
                -9.0168600000e-05, -3.2937400000e-04, -1.0771500000e-03,
                -3.2189200000e-03, -8.8784500000e-03, -2.2554500000e-02,
                -5.1845300000e-02, -1.0372200000e-01, -1.6659500000e-01,
                -1.6016500000e-01, 6.2365400000e-02,  4.6553400000e-01,
                5.8156400000e-01,  1.9136900000e-01,  -1.2856600000e-01,
                -1.8321600000e-01, -5.6626600000e-02},
      doubles_t{7.4010000000e+06, 1.1080000000e+06, 2.5210000000e+05,
                7.1380000000e+04, 2.3260000000e+04, 8.3900000000e+03,
                3.2710000000e+03, 1.3570000000e+03, 5.9200000000e+02,
                2.6910000000e+02, 1.2650000000e+02, 6.1030000000e+01,
                2.9860000000e+01, 1.4170000000e+01, 7.0220000000e+00,
                3.5110000000e+00, 1.7580000000e+00, 7.8410000000e-01,
                3.4800000000e-01, 1.4910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7580000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4800000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.3037900000e-09,  1.7422100000e-08,  9.5530000000e-08,
                3.8361900000e-07,  1.5053500000e-06,  4.5863200000e-06,
                1.5500200000e-05,  4.0669500000e-05,  1.3133700000e-04,
                3.3268500000e-04,  1.0667000000e-03,  2.5721000000e-03,
                4.8358400000e-03,  -3.8572800000e-03, -4.7612500000e-02,
                -1.8295700000e-01, -7.3282000000e-02, 4.3482200000e-01,
                5.8171100000e-01,  1.7937300000e-01},
      doubles_t{7.4010000000e+06, 1.1080000000e+06, 2.5210000000e+05,
                7.1380000000e+04, 2.3260000000e+04, 8.3900000000e+03,
                3.2710000000e+03, 1.3570000000e+03, 5.9200000000e+02,
                2.6910000000e+02, 1.2650000000e+02, 6.1030000000e+01,
                2.9860000000e+01, 1.4170000000e+01, 7.0220000000e+00,
                3.5110000000e+00, 1.7580000000e+00, 7.8410000000e-01,
                3.4800000000e-01, 1.4910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4910000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.8833800000e-04, 1.6428800000e-03, 8.9489300000e-03,
                3.5511100000e-02, 1.0514700000e-01, 2.2155200000e-01,
                3.0841200000e-01, 2.3220300000e-01, 1.2447800000e-01,
                1.6831000000e-01, 1.6829200000e-01, 5.9798100000e-02},
      doubles_t{2.9270000000e+03, 6.9350000000e+02, 2.2470000000e+02,
                8.5170000000e+01, 3.5530000000e+01, 1.5730000000e+01,
                7.1650000000e+00, 3.3220000000e+00, 1.4780000000e+00,
                6.5520000000e-01, 2.7510000000e-01, 1.0970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4780000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5520000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5020700000e-04, -1.3092800000e-03, -7.1650100000e-03,
                -2.8572300000e-02, -8.6015800000e-02, -1.8399600000e-01,
                -2.6584400000e-01, -1.8874100000e-01, 1.2537800000e-01,
                3.9724600000e-01, 3.9465800000e-01, 1.4191700000e-01},
      doubles_t{2.9270000000e+03, 6.9350000000e+02, 2.2470000000e+02,
                8.5170000000e+01, 3.5530000000e+01, 1.5730000000e+01,
                7.1650000000e+00, 3.3220000000e+00, 1.4780000000e+00,
                6.5520000000e-01, 2.7510000000e-01, 1.0970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7510000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0518000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5160000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3320000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6680000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5620000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2640000000e+00}));
    return abs_t(name, 18, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_5_plus_d_cparen_z_18

} // namespace chemcache
