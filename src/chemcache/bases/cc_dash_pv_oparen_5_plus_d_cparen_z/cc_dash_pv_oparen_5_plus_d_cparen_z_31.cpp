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

abs_t cc_dash_pv_oparen_5_plus_d_cparen_z_31() {
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
      doubles_t{2.4000000000e-07,  1.8600000000e-06,  9.8000000000e-06,
                4.1520000000e-05,  1.5205000000e-04,  5.0077000000e-04,
                1.5187000000e-03,  4.3025000000e-03,  1.1452300000e-02,
                2.8564000000e-02,  6.5748500000e-02,  1.3528950000e-01,
                2.3455140000e-01,  3.0783510000e-01,  2.5299470000e-01,
                9.6010400000e-02,  9.7885000000e-03,  5.9120000000e-04,
                -5.5400000000e-05, 1.3800000000e-05,  -6.4200000000e-05,
                1.6900000000e-05,  -1.1600000000e-05, 5.8600000000e-06,
                -2.2900000000e-06, 6.1000000000e-07},
      doubles_t{
        1.0861520000e+08, 1.6264540000e+07, 3.7001120000e+06, 1.0471690000e+06,
        3.4106760000e+05, 1.2277150000e+05, 4.7659580000e+04, 1.9633350000e+04,
        8.4887350000e+03, 3.8231380000e+03, 1.7844760000e+03, 8.6005300000e+02,
        4.2669870000e+02, 2.1726160000e+02, 1.1296990000e+02, 5.9449440000e+01,
        3.0782260000e+01, 1.6423210000e+01, 8.7578890000e+00, 4.4096290000e+00,
        2.2494490000e+00, 1.1261150000e+00, 5.1548600000e-01, 2.4257800000e-01,
        1.0708600000e-01, 4.6988000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.0000000000e-08, -5.8000000000e-07, -3.0300000000e-06,
                -1.2870000000e-05, -4.7140000000e-05, -1.5530000000e-04,
                -4.7180000000e-04, -1.3405000000e-03, -3.5955000000e-03,
                -9.1016000000e-03, -2.1636000000e-02, -4.7336500000e-02,
                -9.2499700000e-02, -1.5043510000e-01, -1.7212270000e-01,
                -4.4017900000e-02, 2.9738280000e-01,  5.2797480000e-01,
                3.0089050000e-01,  4.5881900000e-02,  1.2828000000e-03,
                1.2588000000e-03,  -2.0110000000e-04, 1.1730000000e-04,
                -4.7400000000e-05, 1.1200000000e-05},
      doubles_t{
        1.0861520000e+08, 1.6264540000e+07, 3.7001120000e+06, 1.0471690000e+06,
        3.4106760000e+05, 1.2277150000e+05, 4.7659580000e+04, 1.9633350000e+04,
        8.4887350000e+03, 3.8231380000e+03, 1.7844760000e+03, 8.6005300000e+02,
        4.2669870000e+02, 2.1726160000e+02, 1.1296990000e+02, 5.9449440000e+01,
        3.0782260000e+01, 1.6423210000e+01, 8.7578890000e+00, 4.4096290000e+00,
        2.2494490000e+00, 1.1261150000e+00, 5.1548600000e-01, 2.4257800000e-01,
        1.0708600000e-01, 4.6988000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{3.0000000000e-08,  2.2000000000e-07,  1.1600000000e-06,
                4.9100000000e-06,  1.7980000000e-05,  5.9200000000e-05,
                1.8010000000e-04,  5.1140000000e-04,  1.3740000000e-03,
                3.4818000000e-03,  8.3169000000e-03,  1.8318000000e-02,
                3.6390300000e-02,  6.0808300000e-02,  7.3293900000e-02,
                1.9741600000e-02,  -1.6129700000e-01, -4.0219480000e-01,
                -2.9272480000e-01, 2.7069420000e-01,  6.3597590000e-01,
                3.7024890000e-01,  4.0282500000e-02,  -1.5554000000e-03,
                5.5340000000e-04,  -7.7900000000e-05},
      doubles_t{
        1.0861520000e+08, 1.6264540000e+07, 3.7001120000e+06, 1.0471690000e+06,
        3.4106760000e+05, 1.2277150000e+05, 4.7659580000e+04, 1.9633350000e+04,
        8.4887350000e+03, 3.8231380000e+03, 1.7844760000e+03, 8.6005300000e+02,
        4.2669870000e+02, 2.1726160000e+02, 1.1296990000e+02, 5.9449440000e+01,
        3.0782260000e+01, 1.6423210000e+01, 8.7578890000e+00, 4.4096290000e+00,
        2.2494490000e+00, 1.1261150000e+00, 5.1548600000e-01, 2.4257800000e-01,
        1.0708600000e-01, 4.6988000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1548600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4257800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0708600000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.0000000000e-09, -5.1000000000e-08, -2.7000000000e-07,
                -1.1420000000e-06, -4.1830000000e-06, -1.3781000000e-05,
                -4.1882000000e-05, -1.1902300000e-04, -3.1960000000e-04,
                -8.1070000000e-04, -1.9360000000e-03, -4.2722000000e-03,
                -8.4945000000e-03, -1.4270900000e-02, -1.7268100000e-02,
                -4.7782000000e-03, 3.9492700000e-02,  1.0272000000e-01,
                7.7352900000e-02,  -8.4956500000e-02, -2.2198340000e-01,
                -2.5320890000e-01, -3.2095400000e-02, 4.3100170000e-01,
                5.7776890000e-01,  1.8148420000e-01},
      doubles_t{
        1.0861520000e+08, 1.6264540000e+07, 3.7001120000e+06, 1.0471690000e+06,
        3.4106760000e+05, 1.2277150000e+05, 4.7659580000e+04, 1.9633350000e+04,
        8.4887350000e+03, 3.8231380000e+03, 1.7844760000e+03, 8.6005300000e+02,
        4.2669870000e+02, 2.1726160000e+02, 1.1296990000e+02, 5.9449440000e+01,
        3.0782260000e+01, 1.6423210000e+01, 8.7578890000e+00, 4.4096290000e+00,
        2.2494490000e+00, 1.1261150000e+00, 5.1548600000e-01, 2.4257800000e-01,
        1.0708600000e-01, 4.6988000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6988000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{2.8300000000e-05, 2.5290000000e-04, 1.4686000000e-03,
                           6.5627000000e-03, 2.3802300000e-02, 7.0894500000e-02,
                           1.6763840000e-01, 2.9597540000e-01, 3.4886100000e-01,
                           2.1754960000e-01, 5.2051100000e-02, 3.4378000000e-03,
                           9.8330000000e-04, 3.9100000000e-05, 4.4900000000e-05,
                           -1.9100000000e-05, 5.3000000000e-06},
                 doubles_t{3.2152190000e+04, 7.6093840000e+03, 2.4714740000e+03,
                           9.4606360000e+02, 4.0194710000e+02, 1.8364690000e+02,
                           8.8533260000e+01, 4.4270360000e+01, 2.2723080000e+01,
                           1.1823140000e+01, 6.0421350000e+00, 3.0317540000e+00,
                           1.4933660000e+00, 7.0972700000e-01, 2.4859300000e-01,
                           9.4395000000e-02, 3.5887000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0700000000e-05, -9.5800000000e-05, -5.5820000000e-04,
                -2.5040000000e-03, -9.1996000000e-03, -2.7997300000e-02,
                -6.8874600000e-02, -1.2738430000e-01, -1.5858890000e-01,
                -4.2496800000e-02, 2.4414400000e-01, 4.4591110000e-01,
                3.5295220000e-01, 1.0494460000e-01, 5.2499000000e-03,
                -1.2190000000e-04, 1.2630000000e-04},
      doubles_t{3.2152190000e+04, 7.6093840000e+03, 2.4714740000e+03,
                9.4606360000e+02, 4.0194710000e+02, 1.8364690000e+02,
                8.8533260000e+01, 4.4270360000e+01, 2.2723080000e+01,
                1.1823140000e+01, 6.0421350000e+00, 3.0317540000e+00,
                1.4933660000e+00, 7.0972700000e-01, 2.4859300000e-01,
                9.4395000000e-02, 3.5887000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0972700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4859300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.4395000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.7000000000e-06, 1.5800000000e-05, 9.0800000000e-05,
                4.1200000000e-04, 1.4984000000e-03, 4.6252000000e-03,
                1.1271300000e-02, 2.1321200000e-02, 2.5952300000e-02,
                6.6320000000e-03, -5.0170400000e-02, -8.4297700000e-02,
                -9.0302300000e-02, 1.8472700000e-02, 3.1232250000e-01,
                5.3648160000e-01, 2.8316500000e-01},
      doubles_t{3.2152190000e+04, 7.6093840000e+03, 2.4714740000e+03,
                9.4606360000e+02, 4.0194710000e+02, 1.8364690000e+02,
                8.8533260000e+01, 4.4270360000e+01, 2.2723080000e+01,
                1.1823140000e+01, 6.0421350000e+00, 3.0317540000e+00,
                1.4933660000e+00, 7.0972700000e-01, 2.4859300000e-01,
                9.4395000000e-02, 3.5887000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5887000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{8.9200000000e-05, 8.6250000000e-04, 5.0094000000e-03,
                1.9964900000e-02, 5.8321400000e-02, 1.3168680000e-01,
                2.2186760000e-01, 2.8250590000e-01, 2.8319890000e-01,
                2.1582510000e-01, 1.0524360000e-01, 2.0723200000e-02},
      doubles_t{1.0405050000e+03, 3.1459710000e+02, 1.2278760000e+02,
                5.4760370000e+01, 2.6298940000e+01, 1.3263450000e+01,
                6.8850650000e+00, 3.5795250000e+00, 1.8315640000e+00,
                9.1290900000e-01, 4.3534000000e-01, 1.8851800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1290900000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3534000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8851800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.5800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8260000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1460000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9860000000e-01}));
    return abs_t(name, 31, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_5_plus_d_cparen_z_31

} // namespace chemcache
