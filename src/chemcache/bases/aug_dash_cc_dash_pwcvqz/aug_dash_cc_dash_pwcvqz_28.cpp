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

#include "aug_dash_cc_dash_pwcvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pwcvqz_28() {
    // Basis Set name and origin point
    std::string name("aug-cc-pwcvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{
        2.0850230000e-06,  1.6214250000e-05, 8.5307730000e-05, 3.6037930000e-04,
        1.3126320000e-03,  4.2723950000e-03, 1.2659940000e-02, 3.4230650000e-02,
        8.3267350000e-02,  1.7457580000e-01, 2.8980230000e-01, 3.2073230000e-01,
        1.7782650000e-01,  3.6330070000e-02, 2.0894580000e-02, 1.5969540000e-02,
        1.7302450000e-03,  4.8173650000e-04, 1.2468000000e-04, 2.3403780000e-05,
        -1.3489620000e-05, 4.7995390000e-06},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.2981550000e-07, -5.6777680000e-06, -2.9859220000e-05,
                -1.2627840000e-04, -4.6002520000e-04, -1.5032300000e-03,
                -4.4791710000e-03, -1.2325120000e-02, -3.0976390000e-02,
                -6.9986510000e-02, -1.3392560000e-01, -1.9668810000e-01,
                -1.4022580000e-01, 1.9607300000e-01,  5.5645170000e-01,
                3.7021670000e-01,  5.3960340000e-02,  4.7064920000e-03,
                3.7361150000e-03,  -7.4680500000e-05, 9.5753360000e-05,
                -3.4371150000e-05},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4532370000e-07,  1.9096370000e-06,  1.0034890000e-05,
                4.2488350000e-05,  1.5461870000e-04,  5.0645500000e-04,
                1.5088400000e-03,  4.1788160000e-03,  1.0566680000e-02,
                2.4393000000e-02,  4.8121010000e-02,  7.5651700000e-02,
                5.7485340000e-02,  -9.7831860000e-02, -3.8645360000e-01,
                -3.6065980000e-01, 2.8886500000e-01,  6.8309170000e-01,
                3.2067780000e-01,  1.6514710000e-02,  3.1246930000e-04,
                5.3476750000e-03},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.1142610000e-08, -3.9807790000e-07, -2.0920330000e-06,
                -8.8568360000e-06, -3.2236020000e-05, -1.0557570000e-04,
                -3.1463560000e-04, -8.7141130000e-04, -2.2054150000e-03,
                -5.0955170000e-03, -1.0082370000e-02, -1.5920700000e-02,
                -1.2220370000e-02, 2.1158690000e-02,  8.7473040000e-02,
                8.5116050000e-02,  -8.1156830000e-02, -2.1663590000e-01,
                -2.0769580000e-01, 2.6992960000e-01,  5.8427110000e-01,
                2.9579030000e-01},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7717330000e-07, -1.3976150000e-06, -7.2038590000e-06,
                -3.1285520000e-05, -1.1026990000e-04, -3.7557350000e-04,
                -1.0679000000e-03, -3.1284270000e-03, -7.4102720000e-03,
                -1.8595010000e-02, -3.3227980000e-02, -6.1184400000e-02,
                -3.2863340000e-02, 4.8442020000e-02,  3.9898800000e-01,
                3.3792900000e-01,  -7.0162240000e-01, -1.6944340000e+00,
                2.5826800000e+00,  -5.4987580000e-01, -1.3945510000e+00,
                1.3616180000e+00},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.7620690000e-08, -7.7411420000e-07, -3.9597460000e-06,
                -1.7369140000e-05, -6.0449530000e-05, -2.0905550000e-04,
                -5.8343120000e-04, -1.7467130000e-03, -4.0267840000e-03,
                -1.0423620000e-02, -1.7811720000e-02, -3.4669960000e-02,
                -1.5409450000e-02, 2.0864550000e-02,  2.2068200000e-01,
                1.0739660000e-01,  -8.3106680000e-02, -8.9429990000e-01,
                3.7501410000e-01,  1.4967530000e+00,  -5.3203900000e-01,
                -8.2639390000e-01},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3701000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.9883110000e-07, -1.4497160000e-06, -8.3634440000e-06,
                -3.1255320000e-05, -1.3276240000e-04, -3.5896620000e-04,
                -1.3414880000e-03, -2.8240180000e-03, -9.8636070000e-03,
                -1.5234810000e-02, -4.9844890000e-02, -3.5223540000e-02,
                -1.0564240000e-01, 2.1380340000e-01,  1.5176750000e-01,
                1.1804280000e+00,  -3.5264900000e+00, 2.5508150000e+00,
                4.3064610000e-01,  -3.4500410000e+00, 4.7403420000e+00,
                -2.3121560000e+00},
      doubles_t{
        1.5503020000e+07, 2.3210710000e+06, 5.2818020000e+05, 1.4960400000e+05,
        4.8809070000e+04, 1.7622390000e+04, 6.8742010000e+03, 2.8518570000e+03,
        1.2443340000e+03, 5.6615470000e+02, 2.6678730000e+02, 1.2937520000e+02,
        6.3651610000e+01, 2.9333020000e+01, 1.4822120000e+01, 7.4605130000e+00,
        3.4135760000e+00, 1.6328840000e+00, 7.4357900000e-01, 1.7259200000e-01,
        7.7420000000e-02, 3.3701000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5001000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7599000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4670000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.1000000000e-05, 9.6000000000e-05, 5.6400000000e-04,
                2.5700000000e-03, 9.6770000000e-03, 3.0742000000e-02,
                8.2065000000e-02, 1.7684500000e-01, 2.9116200000e-01,
                3.3081600000e-01, 2.0881900000e-01, 5.3739000000e-02,
                1.7090000000e-03, -1.7980000000e-03, -1.0460000000e-03,
                -6.2000000000e-05, -1.0000000000e-06},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 3.5643000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-4.0000000000e-06, -3.5000000000e-05, -2.0200000000e-04,
                -9.2500000000e-04, -3.4990000000e-03, -1.1280000000e-02,
                -3.0819000000e-02, -6.9036000000e-02, -1.1873300000e-01,
                -1.4244300000e-01, -4.3107000000e-02, 2.0349900000e-01,
                4.0622800000e-01, 3.7093800000e-01, 1.6431200000e-01,
                1.7265000000e-02, -1.6820000000e-03, 6.5100000000e-04},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 9.1546000000e-02, 3.5643000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06, 1.4000000000e-05, 8.4000000000e-05,
                3.7100000000e-04, 1.4600000000e-03, 4.5370000000e-03,
                1.2922000000e-02, 2.7766000000e-02, 5.0956000000e-02,
                5.6745000000e-02, 2.4119000000e-02, -1.2508700000e-01,
                -2.0654100000e-01, -2.1498300000e-01, 2.0781900000e-01,
                6.3854200000e-01, 3.2628400000e-01, 1.2981000000e-02},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 9.1546000000e-02, 3.5643000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0000000000e-06, 3.6000000000e-05, 1.9900000000e-04,
                9.5900000000e-04, 3.4630000000e-03, 1.1756000000e-02,
                3.0690000000e-02, 7.3081000000e-02, 1.1682700000e-01,
                1.9116300000e-01, 1.0654300000e-01, -9.4992400000e-01,
                -5.6360500000e-01, 1.9346110000e+00, -5.1017700000e-01,
                -1.1722970000e+00, 1.0158980000e+00, 2.3798000000e-02},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 9.1546000000e-02, 3.5643000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06, 2.4000000000e-05, 1.5500000000e-04,
                6.4600000000e-04, 2.6890000000e-03, 7.8990000000e-03,
                2.3912000000e-02, 4.8006000000e-02, 9.6213000000e-02,
                9.7582000000e-02, 7.7022000000e-02, -3.6013600000e-01,
                -4.4558900000e-01, -1.1653800000e-01, 1.1705030000e+00,
                -3.8931000000e-02, -7.6881400000e-01, -5.9410000000e-03},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 9.1546000000e-02, 3.5643000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 7.0000000000e-06, 4.1000000000e-05,
                1.8600000000e-04, 7.1200000000e-04, 2.2690000000e-03,
                6.2910000000e-03, 1.3934000000e-02, 2.4531000000e-02,
                2.8770000000e-02, 8.7810000000e-03, -5.1731000000e-02,
                -9.5945000000e-02, -1.1067500000e-01, 2.7270000000e-03,
                3.1502700000e-01, 5.4281000000e-01, 2.7575200000e-01},
      doubles_t{4.5066620000e+04, 1.0662640000e+04, 3.4634310000e+03,
                1.3266240000e+03, 5.6474140000e+02, 2.5882880000e+02,
                1.2531150000e+02, 6.3242330000e+01, 3.2824060000e+01,
                1.7396920000e+01, 9.3243940000e+00, 4.9388390000e+00,
                2.5620250000e+00, 1.2992940000e+00, 6.3467800000e-01,
                2.3226200000e-01, 9.1546000000e-02, 3.5643000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5643000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0416500000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9093000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3880000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9565000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6674000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.5000000000e-04, 3.2000000000e-03, 1.6477000000e-02,
                           5.6001000000e-02, 1.3923400000e-01, 2.4277600000e-01,
                           3.0406900000e-01, 2.9331900000e-01, 2.1725900000e-01,
                           1.0469900000e-01, 1.7339000000e-02},
                 doubles_t{4.1434270000e+02, 1.2439260000e+02, 4.8115400000e+01,
                           2.0842000000e+01, 9.6221000000e+00, 4.6034000000e+00,
                           2.1963000000e+00, 1.0206000000e+00, 4.5440000000e-01,
                           1.9000000000e-01, 7.1500000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.6200000000e-04, -3.3350000000e-03, -1.7176000000e-02,
                -5.9245000000e-02, -1.5123200000e-01, -2.5446100000e-01,
                -2.2871800000e-01, 2.0616000000e-02, 3.4275300000e-01,
                4.7387100000e-01, 1.8343500000e-01},
      doubles_t{4.1434270000e+02, 1.2439260000e+02, 4.8115400000e+01,
                2.0842000000e+01, 9.6221000000e+00, 4.6034000000e+00,
                2.1963000000e+00, 1.0206000000e+00, 4.5440000000e-01,
                1.9000000000e-01, 7.1500000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-5.2900000000e-04, -4.9280000000e-03, -2.5372000000e-02,
                -8.9654000000e-02, -2.3874400000e-01, -3.5972100000e-01,
                4.4440000000e-03, 5.9654700000e-01, 3.3130200000e-01,
                -5.0900900000e-01, -3.5013900000e-01},
      doubles_t{4.1434270000e+02, 1.2439260000e+02, 4.8115400000e+01,
                2.0842000000e+01, 9.6221000000e+00, 4.6034000000e+00,
                2.1963000000e+00, 1.0206000000e+00, 4.5440000000e-01,
                1.9000000000e-01, 7.1500000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-7.8700000000e-04, -6.7780000000e-03, -3.7726000000e-02,
                -1.2912500000e-01, -4.1480900000e-01, -2.8395700000e-01,
                8.6719400000e-01, 2.7232700000e-01, -9.8684200000e-01,
                1.8536000000e-01, 5.0261000000e-01},
      doubles_t{4.1434270000e+02, 1.2439260000e+02, 4.8115400000e+01,
                2.0842000000e+01, 9.6221000000e+00, 4.6034000000e+00,
                2.1963000000e+00, 1.0206000000e+00, 4.5440000000e-01,
                1.9000000000e-01, 7.1500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6910000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0379100000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4257000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8872000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5380000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2650000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5105000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8307000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5457000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7448000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3495000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4486000000e+00}));
    return abs_t(name, 28, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pwcvqz_28

} // namespace chemcache
