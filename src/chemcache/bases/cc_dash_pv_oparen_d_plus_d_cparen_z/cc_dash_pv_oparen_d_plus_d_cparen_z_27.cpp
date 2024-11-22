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

#include "cc_dash_pv_oparen_d_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_d_plus_d_cparen_z_27() {
    // Basis Set name and origin point
    std::string name("cc-pv(d+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.9790260000e-06,  6.2040710000e-05, 3.2617350000e-04,
                1.3753600000e-03,  4.9799970000e-03, 1.5964340000e-02,
                4.5520860000e-02,  1.1273850000e-01, 2.2682620000e-01,
                3.2030740000e-01,  2.3740210000e-01, 7.4776860000e-02,
                9.5818720000e-02,  9.6499110000e-02, 1.6233620000e-02,
                -4.5354970000e-04, 5.1135190000e-05, -4.1745080000e-05,
                4.0275770000e-05,  -5.7890670000e-06},
      doubles_t{4.6756750000e+06, 7.0016150000e+05, 1.5933730000e+05,
                4.5130460000e+04, 1.4722380000e+04, 5.3142220000e+03,
                2.0720180000e+03, 8.5861880000e+02, 3.7354970000e+02,
                1.6892290000e+02, 7.8296390000e+01, 3.5521230000e+01,
                1.7041440000e+01, 8.1730000000e+00, 3.6103180000e+00,
                1.6970470000e+00, 7.4353200000e-01, 1.5834400000e-01,
                7.5036000000e-02, 3.3091000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.2002400000e-06, -3.2658310000e-05, -1.7176440000e-04,
                -7.2478530000e-04, -2.6314620000e-03, -8.4892720000e-03,
                -2.4606190000e-02, -6.3220590000e-02, -1.3819570000e-01,
                -2.3406800000e-01, -2.4150020000e-01, 3.0353120000e-02,
                5.1013410000e-01,  4.9749390000e-01,  8.9707460000e-02,
                -5.9410340000e-03, 2.1753620000e-04,  -5.4801550000e-04,
                4.5258040000e-04,  -1.0667480000e-04},
      doubles_t{4.6756750000e+06, 7.0016150000e+05, 1.5933730000e+05,
                4.5130460000e+04, 1.4722380000e+04, 5.3142220000e+03,
                2.0720180000e+03, 8.5861880000e+02, 3.7354970000e+02,
                1.6892290000e+02, 7.8296390000e+01, 3.5521230000e+01,
                1.7041440000e+01, 8.1730000000e+00, 3.6103180000e+00,
                1.6970470000e+00, 7.4353200000e-01, 1.5834400000e-01,
                7.5036000000e-02, 3.3091000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.5926920000e-07,  7.4618510000e-06,  3.9231370000e-05,
                1.6577060000e-04,  6.0243350000e-04,  1.9552170000e-03,
                5.7263260000e-03,  1.5129840000e-02,  3.4839730000e-02,
                6.5703510000e-02,  7.8315030000e-02,  -1.8770370000e-02,
                -3.0626630000e-01, -4.5664290000e-01, 1.3781690000e-01,
                7.1936760000e-01,  3.9925790000e-01,  2.0799330000e-02,
                -7.8206630000e-03, 3.5339110000e-03},
      doubles_t{4.6756750000e+06, 7.0016150000e+05, 1.5933730000e+05,
                4.5130460000e+04, 1.4722380000e+04, 5.3142220000e+03,
                2.0720180000e+03, 8.5861880000e+02, 3.7354970000e+02,
                1.6892290000e+02, 7.8296390000e+01, 3.5521230000e+01,
                1.7041440000e+01, 8.1730000000e+00, 3.6103180000e+00,
                1.6970470000e+00, 7.4353200000e-01, 1.5834400000e-01,
                7.5036000000e-02, 3.3091000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.0288400000e-07, -1.5775800000e-06, -8.2988130000e-06,
                -3.5041540000e-05, -1.2746550000e-04, -4.1326950000e-04,
                -1.2122610000e-03, -3.1993180000e-03, -7.3909720000e-03,
                -1.3936490000e-02, -1.6785750000e-02, 4.1498560000e-03,
                6.7976460000e-02,  1.0758070000e-01,  -4.1660220000e-02,
                -2.1280440000e-01, -2.3813600000e-01, 2.6507880000e-01,
                5.7227740000e-01,  3.0915560000e-01},
      doubles_t{4.6756750000e+06, 7.0016150000e+05, 1.5933730000e+05,
                4.5130460000e+04, 1.4722380000e+04, 5.3142220000e+03,
                2.0720180000e+03, 8.5861880000e+02, 3.7354970000e+02,
                1.6892290000e+02, 7.8296390000e+01, 3.5521230000e+01,
                1.7041440000e+01, 8.1730000000e+00, 3.6103180000e+00,
                1.6970470000e+00, 7.4353200000e-01, 1.5834400000e-01,
                7.5036000000e-02, 3.3091000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.8630530000e-07, -3.0687880000e-06, -1.5648260000e-05,
                -6.8835880000e-05, -2.3773670000e-04, -8.2131730000e-04,
                -2.2296300000e-03, -6.4678410000e-03, -1.3254630000e-02,
                -2.9466860000e-02, -2.5990660000e-02, -8.4998070000e-03,
                1.7273160000e-01,  1.5121890000e-01,  3.5545090000e-02,
                -8.8293530000e-01, 2.1435300000e-01,  1.7118650000e+00,
                -7.1400370000e-01, -8.0277270000e-01},
      doubles_t{4.6756750000e+06, 7.0016150000e+05, 1.5933730000e+05,
                4.5130460000e+04, 1.4722380000e+04, 5.3142220000e+03,
                2.0720180000e+03, 8.5861880000e+02, 3.7354970000e+02,
                1.6892290000e+02, 7.8296390000e+01, 3.5521230000e+01,
                1.7041440000e+01, 8.1730000000e+00, 3.6103180000e+00,
                1.6970470000e+00, 7.4353200000e-01, 1.5834400000e-01,
                7.5036000000e-02, 3.3091000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3091000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.1000000000e-05, 3.6900000000e-04, 2.1280000000e-03,
                9.3720000000e-03, 3.3155000000e-02, 9.4752000000e-02,
                2.0909300000e-01, 3.3372200000e-01, 3.3220800000e-01,
                1.5461300000e-01, 2.0902000000e-02, -2.0240000000e-03,
                -1.6970000000e-03, -2.8000000000e-04, 2.6000000000e-05,
                -1.0000000000e-05},
      doubles_t{1.9267780000e+04, 4.5609860000e+03, 1.4814360000e+03,
                5.6686710000e+02, 2.4049100000e+02, 1.0961050000e+02,
                5.2594910000e+01, 2.6083610000e+01, 1.3261430000e+01,
                6.7997780000e+00, 3.3934140000e+00, 1.6487660000e+00,
                7.7628200000e-01, 2.9800300000e-01, 1.1361800000e-01,
                4.1624000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5000000000e-05, -1.3100000000e-04, -7.5800000000e-04,
                -3.3630000000e-03, -1.2054000000e-02, -3.5424000000e-02,
                -8.1287000000e-02, -1.3690800000e-01, -1.3901900000e-01,
                3.5468000000e-02, 3.3849800000e-01, 4.5443300000e-01,
                2.7979300000e-01, 4.4776000000e-02, -3.1510000000e-03,
                1.3170000000e-03},
      doubles_t{1.9267780000e+04, 4.5609860000e+03, 1.4814360000e+03,
                5.6686710000e+02, 2.4049100000e+02, 1.0961050000e+02,
                5.2594910000e+01, 2.6083610000e+01, 1.3261430000e+01,
                6.7997780000e+00, 3.3934140000e+00, 1.6487660000e+00,
                7.7628200000e-01, 2.9800300000e-01, 1.1361800000e-01,
                4.1624000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.0000000000e-06, 4.5000000000e-05, 2.5500000000e-04,
                1.1440000000e-03, 4.0610000000e-03, 1.2095000000e-02,
                2.7476000000e-02, 4.7557000000e-02, 4.7302000000e-02,
                -1.4418000000e-02, -1.5006200000e-01, -1.9909200000e-01,
                -7.9783000000e-02, 4.5903500000e-01, 6.1749500000e-01,
                6.4690000000e-02},
      doubles_t{1.9267780000e+04, 4.5609860000e+03, 1.4814360000e+03,
                5.6686710000e+02, 2.4049100000e+02, 1.0961050000e+02,
                5.2594910000e+01, 2.6083610000e+01, 1.3261430000e+01,
                6.7997780000e+00, 3.3934140000e+00, 1.6487660000e+00,
                7.7628200000e-01, 2.9800300000e-01, 1.1361800000e-01,
                4.1624000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.0000000000e-06, -2.9000000000e-05, -1.6700000000e-04,
                -7.4200000000e-04, -2.6620000000e-03, -7.8410000000e-03,
                -1.8051000000e-02, -3.0580000000e-02, -3.1312000000e-02,
                1.1311000000e-02, 8.9990000000e-02, 1.3073300000e-01,
                7.1808000000e-02, -2.2165800000e-01, -5.7102500000e-01,
                -3.6378900000e-01},
      doubles_t{1.9267780000e+04, 4.5609860000e+03, 1.4814360000e+03,
                5.6686710000e+02, 2.4049100000e+02, 1.0961050000e+02,
                5.2594910000e+01, 2.6083610000e+01, 1.3261430000e+01,
                6.7997780000e+00, 3.3934140000e+00, 1.6487660000e+00,
                7.7628200000e-01, 2.9800300000e-01, 1.1361800000e-01,
                4.1624000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1624000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.5100000000e-03, 2.5884000000e-02, 1.0005800000e-01,
                           2.4054700000e-01, 3.5684300000e-01, 3.5957900000e-01,
                           2.3662900000e-01, 6.2129000000e-02},
                 doubles_t{1.2626400000e+02, 3.7522600000e+01, 1.3802100000e+01,
                           5.6092700000e+00, 2.3336900000e+00, 9.3641500000e-01,
                           3.4823700000e-01, 1.1235300000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-4.0670000000e-03, -3.0053000000e-02, -1.1962000000e-01,
                -2.9151300000e-01, -3.1804800000e-01, 9.1698000000e-02,
                5.6082300000e-01, 3.5867800000e-01},
      doubles_t{1.2626400000e+02, 3.7522600000e+01, 1.3802100000e+01,
                5.6092700000e+00, 2.3336900000e+00, 9.3641500000e-01,
                3.4823700000e-01, 1.1235300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1235300000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3, doubles_t{4.2396600000e-01, 7.6842900000e-01},
                 doubles_t{3.7724000000e+00, 9.1700000000e-01}));
    return abs_t(name, 27, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_d_plus_d_cparen_z_27

} // namespace chemcache
