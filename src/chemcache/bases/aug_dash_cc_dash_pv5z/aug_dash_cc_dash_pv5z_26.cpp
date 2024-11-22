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

#include "aug_dash_cc_dash_pv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv5z_26() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{8.6575600000e-08,  5.7791200000e-07, 3.3516450000e-06,
                1.5357650000e-05,  5.9496900000e-05, 2.0141520000e-04,
                6.1490370000e-04,  1.7313210000e-03, 4.5693680000e-03,
                1.1397910000e-02,  2.6867840000e-02, 5.9107620000e-02,
                1.1784260000e-01,  2.0137620000e-01, 2.6607740000e-01,
                2.2598620000e-01,  9.8938900000e-02, 4.9325590000e-02,
                8.6043760000e-02,  7.3224620000e-02, 2.1174500000e-02,
                1.6231430000e-03,  1.1972420000e-06, -5.9636730000e-05,
                -7.9818020000e-07, 9.6647110000e-06, 9.1322000000e-06,
                3.5540720000e-06},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.4668330000e-08, -2.9817080000e-07, -1.7292840000e-06,
                -7.9239330000e-06, -3.0700410000e-05, -1.0394950000e-04,
                -3.1750450000e-04, -8.9494870000e-04, -2.3678120000e-03,
                -5.9369260000e-03, -1.4144860000e-02, -3.1784130000e-02,
                -6.6046300000e-02, -1.2227030000e-01, -1.8869200000e-01,
                -2.0909700000e-01, -9.6988960000e-02, 1.7319160000e-01,
                4.6247720000e-01,  3.9983400000e-01,  1.1524080000e-01,
                8.4087470000e-03,  -1.2969100000e-03, -7.5502380000e-04,
                -6.9074160000e-05, 8.7776040000e-05,  3.2507020000e-05,
                2.2254470000e-05},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0254710000e-08,  6.8466790000e-08,  3.9697900000e-07,
                1.8196510000e-06,  7.0476510000e-06,  2.3877740000e-05,
                7.2914720000e-05,  2.0583920000e-04,  5.4510350000e-04,
                1.3732780000e-03,  3.2931170000e-03,  7.5207460000e-03,
                1.6046520000e-02,  3.1302650000e-02,  5.2528650000e-02,
                6.6052950000e-02,  3.2493770000e-02,  -8.4768090000e-02,
                -2.9501800000e-01, -3.7095470000e-01, -3.4648130000e-02,
                4.4624860000e-01,  5.7175860000e-01,  2.3486780000e-01,
                1.5062330000e-02,  -3.4373390000e-03, 3.2541880000e-03,
                1.3460400000e-04},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.2174270000e-09, -1.4805950000e-08, -8.5839010000e-08,
                -3.9350700000e-07, -1.5238830000e-06, -5.1638350000e-06,
                -1.5765580000e-05, -4.4517920000e-05, -1.1786200000e-04,
                -2.9705060000e-04, -7.1215000000e-04, -1.6275780000e-03,
                -3.4728130000e-03, -6.7863870000e-03, -1.1404470000e-02,
                -1.4420200000e-02, -7.1121990000e-03, 1.8751960000e-02,
                6.8000670000e-02,  8.8975490000e-02,  5.3756650000e-03,
                -1.2624550000e-01, -2.0477760000e-01, -1.8133730000e-01,
                8.9267550000e-02,  4.1523510000e-01,  4.9804810000e-01,
                1.5983440000e-01},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.4248240000e-09, -2.9307920000e-08, -1.7167100000e-07,
                -7.7733200000e-07, -3.0560400000e-06, -1.0164530000e-05,
                -3.1751120000e-05, -8.7198080000e-05, -2.3862280000e-04,
                -5.7862080000e-04, -1.4502430000e-03, -3.1508620000e-03,
                -7.1326340000e-03, -1.2990590000e-02, -2.3997570000e-02,
                -2.6064630000e-02, -2.0372160000e-02, 5.0156810000e-02,
                1.1719660000e-01,  2.3763840000e-01,  -7.2873530000e-02,
                -2.4958890000e-01, -7.3172220000e-01, 2.5624260000e-01,
                1.3560610000e+00,  5.8740610000e-02,  -5.6239170000e-01,
                -5.7462050000e-01},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.3106110000e-09, -4.8241910000e-08, -2.8392120000e-07,
                -1.2782900000e-06, -5.0606590000e-06, -1.6687500000e-05,
                -5.2680540000e-05, -1.4282690000e-04, -3.9688420000e-04,
                -9.4534970000e-04, -2.4190170000e-03, -5.1349810000e-03,
                -1.1952440000e-02, -2.1093590000e-02, -4.0758240000e-02,
                -4.1371040000e-02, -3.8733700000e-02, 9.3455790000e-02,
                1.9220050000e-01,  4.8980340000e-01,  -2.4208920000e-01,
                -8.7646170000e-01, -1.0744520000e+00, 2.3517730000e+00,
                5.2404340000e-01,  -2.2517260000e+00, 2.6958950000e-01,
                9.3682970000e-01},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0545760000e-08, -7.0562980000e-08, -4.0797800000e-07,
                -1.8765380000e-06, -7.2370940000e-06, -2.4649740000e-05,
                -7.4789650000e-05, -2.1282810000e-04, -5.5852570000e-04,
                -1.4234200000e-03, -3.3754760000e-03, -7.8381240000e-03,
                -1.6541570000e-02, -3.3193640000e-02, -5.5331500000e-02,
                -7.4317180000e-02, -3.5762610000e-02, 9.5552640000e-02,
                4.9204630000e-01,  1.2092390000e+00,  -2.5331060000e+00,
                -2.3277190000e+00, 8.1954170000e+00,  -7.2448810000e+00,
                3.1292460000e+00,  1.1613620000e+00,  -3.8758240000e+00,
                2.0596130000e+00},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.0150170000e-09,  6.1066750000e-08,  3.4758290000e-07,
                1.6288970000e-06,  6.1397420000e-06,  2.1508510000e-05,
                6.3029670000e-05,  1.8703180000e-04,  4.6676230000e-04,
                1.2603980000e-03,  2.7936580000e-03,  6.9929300000e-03,
                1.3473230000e-02,  2.9962880000e-02,  4.3040640000e-02,
                7.1140080000e-02,  9.9112400000e-03,  -3.6357680000e-02,
                -4.4357270000e-01, -5.1596470000e-01, 3.2765210000e-01,
                2.9003630000e+00,  -2.8382230000e+00, -1.0430590000e+00,
                4.5316650000e+00,  -3.2699620000e+00, -1.0260610000e+00,
                1.4622030000e+00},
      doubles_t{1.4976900000e+08, 2.6722500000e+07, 5.9590000000e+06,
                1.5936100000e+06, 4.9346400000e+05, 1.7179000000e+05,
                6.5620900000e+04, 2.6965600000e+04, 1.1735300000e+04,
                5.3434300000e+03, 2.5225000000e+03, 1.2264300000e+03,
                6.1126900000e+02, 3.1130800000e+02, 1.6162400000e+02,
                8.5376700000e+01, 4.5791100000e+01, 2.4867700000e+01,
                1.3620200000e+01, 7.4810600000e+00, 4.0881700000e+00,
                2.1990500000e+00, 1.1481500000e+00, 5.7726500000e-01,
                1.9763500000e-01, 1.0898600000e-01, 5.3216000000e-02,
                2.5571000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5571000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2290000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06,  2.2000000000e-05,  1.3700000000e-04,
                6.3300000000e-04,  2.4280000000e-03,  8.0880000000e-03,
                2.3801000000e-02,  6.1490000000e-02,  1.3389600000e-01,
                2.3159000000e-01,  3.0516400000e-01,  2.7096100000e-01,
                1.2908500000e-01,  2.3432000000e-02,  -9.1200000000e-04,
                -1.8800000000e-03, -7.7400000000e-04, -5.0000000000e-05,
                -2.0000000000e-06, -1.0000000000e-06},
      doubles_t{8.6947300000e+04, 2.1254000000e+04, 6.7411400000e+03,
                2.5571600000e+03, 1.0935100000e+03, 5.0616100000e+02,
                2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
                3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01,
                6.0435900000e+00, 3.3024500000e+00, 1.7638900000e+00,
                9.2138600000e-01, 4.6037300000e-01, 1.8361300000e-01,
                7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0000000000e-06, -8.0000000000e-06, -4.8000000000e-05,
                -2.2300000000e-04, -8.5800000000e-04, -2.8740000000e-03,
                -8.5490000000e-03, -2.2528000000e-02, -5.0612000000e-02,
                -9.0853000000e-02, -1.2616700000e-01, -1.0878300000e-01,
                3.8014000000e-02,  2.6860100000e-01,  4.0111800000e-01,
                3.2015900000e-01,  1.1782800000e-01,  1.0321000000e-02,
                -2.0700000000e-04, 2.1000000000e-04},
      doubles_t{8.6947300000e+04, 2.1254000000e+04, 6.7411400000e+03,
                2.5571600000e+03, 1.0935100000e+03, 5.0616100000e+02,
                2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
                3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01,
                6.0435900000e+00, 3.3024500000e+00, 1.7638900000e+00,
                9.2138600000e-01, 4.6037300000e-01, 1.8361300000e-01,
                7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  1.0000000000e-05,  5.7000000000e-05,
                2.9800000000e-04,  1.0240000000e-03,  3.8480000000e-03,
                1.0319000000e-02,  3.0711000000e-02,  6.1799000000e-02,
                1.3000200000e-01,  1.5210100000e-01,  3.8027400000e-01,
                -8.1249700000e-01, -1.5336760000e+00, 2.7453080000e+00,
                -1.9352200000e-01, -2.1517710000e+00, 2.2005390000e+00,
                -1.0191770000e+00, -9.2503000000e-02},
      doubles_t{8.6947300000e+04, 2.1254000000e+04, 6.7411400000e+03,
                2.5571600000e+03, 1.0935100000e+03, 5.0616100000e+02,
                2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
                3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01,
                6.0435900000e+00, 3.3024500000e+00, 1.7638900000e+00,
                9.2138600000e-01, 4.6037300000e-01, 1.8361300000e-01,
                7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        3.0000000000e-06, 1.6000000000e-05, 7.2000000000e-05, 2.8500000000e-04,
        9.3100000000e-04, 2.8370000000e-03, 7.3090000000e-03, 1.6890000000e-02,
        2.9488000000e-02, 4.3315000000e-02, 3.4262000000e-02, -1.2642000000e-02,
        -1.1538000000e-01, -1.5979900000e-01, -1.6246600000e-01,
        1.1574900000e-01, 6.3340800000e-01, 3.8482200000e-01, 1.3843000000e-02},
      doubles_t{
        2.1254000000e+04, 6.7411400000e+03, 2.5571600000e+03, 1.0935100000e+03,
        5.0616100000e+02, 2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
        3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01, 6.0435900000e+00,
        3.3024500000e+00, 1.7638900000e+00, 9.2138600000e-01, 4.6037300000e-01,
        1.8361300000e-01, 7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  5.0000000000e-06,  3.5000000000e-05,
                1.4800000000e-04,  6.3000000000e-04,  1.9130000000e-03,
                6.2820000000e-03,  1.5025000000e-02,  3.7670000000e-02,
                5.9521000000e-02,  1.0230200000e-01,  6.6716000000e-02,
                -4.6980000000e-03, -3.9425400000e-01, -4.1024800000e-01,
                -1.0141000000e-02, 1.1891730000e+00,  -4.7688000000e-02,
                -7.4062300000e-01, -1.3427000000e-02},
      doubles_t{8.6947300000e+04, 2.1254000000e+04, 6.7411400000e+03,
                2.5571600000e+03, 1.0935100000e+03, 5.0616100000e+02,
                2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
                3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01,
                6.0435900000e+00, 3.3024500000e+00, 1.7638900000e+00,
                9.2138600000e-01, 4.6037300000e-01, 1.8361300000e-01,
                7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  8.0000000000e-06,  3.9000000000e-05,
                2.2400000000e-04,  6.9300000000e-04,  2.8800000000e-03,
                6.9440000000e-03,  2.2747000000e-02,  4.0625000000e-02,
                9.6658000000e-02,  9.2204000000e-02,  2.1329700000e-01,
                -2.0605400000e-01, -6.9704000000e-01, -5.8617900000e-01,
                2.0953970000e+00,  -6.2060400000e-01, -1.2761030000e+00,
                1.0949440000e+00,  3.0999000000e-02},
      doubles_t{8.6947300000e+04, 2.1254000000e+04, 6.7411400000e+03,
                2.5571600000e+03, 1.0935100000e+03, 5.0616100000e+02,
                2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
                3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01,
                6.0435900000e+00, 3.3024500000e+00, 1.7638900000e+00,
                9.2138600000e-01, 4.6037300000e-01, 1.8361300000e-01,
                7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        2.0000000000e-06, 1.1000000000e-05, 4.9000000000e-05, 1.9000000000e-04,
        6.2900000000e-04, 1.8960000000e-03, 4.9400000000e-03, 1.1284000000e-02,
        1.9975000000e-02, 2.8652000000e-02, 2.3447000000e-02, -1.0157000000e-02,
        -7.3781000000e-02, -1.0603300000e-01, -1.0721400000e-01,
        3.8887000000e-02, 3.5600200000e-01, 5.2561900000e-01, 2.2889900000e-01},
      doubles_t{
        2.1254000000e+04, 6.7411400000e+03, 2.5571600000e+03, 1.0935100000e+03,
        5.0616100000e+02, 2.4713900000e+02, 1.2541500000e+02, 6.5651200000e+01,
        3.5307900000e+01, 1.9426000000e+01, 1.0843800000e+01, 6.0435900000e+00,
        3.3024500000e+00, 1.7638900000e+00, 9.2138600000e-01, 4.6037300000e-01,
        1.8361300000e-01, 7.5857000000e-02, 3.1147000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1147000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2790000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.6100000000e-04, 1.5040000000e-03, 8.2640000000e-03,
                3.0308000000e-02, 8.3310000000e-02, 1.7347700000e-01,
                2.6148300000e-01, 2.9968000000e-01, 2.6956400000e-01,
                1.7862600000e-01, 6.9195000000e-02, 7.4640000000e-03},
      doubles_t{4.8911500000e+02, 1.4720800000e+02, 5.7127800000e+01,
                2.5042900000e+01, 1.1667200000e+01, 5.6825000000e+00,
                2.8132500000e+00, 1.3773300000e+00, 6.5653400000e-01,
                3.0070000000e-01, 1.3018700000e-01, 5.1582000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-1.7700000000e-04, -1.6620000000e-03, -9.1670000000e-03,
                -3.3813000000e-02, -9.5253000000e-02, -2.0097300000e-01,
                -2.7166500000e-01, -1.6189900000e-01, 1.6333400000e-01,
                4.3601600000e-01, 3.6887400000e-01, 7.7226000000e-02},
      doubles_t{4.8911500000e+02, 1.4720800000e+02, 5.7127800000e+01,
                2.5042900000e+01, 1.1667200000e+01, 5.6825000000e+00,
                2.8132500000e+00, 1.3773300000e+00, 6.5653400000e-01,
                3.0070000000e-01, 1.3018700000e-01, 5.1582000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.5500000000e-04, 2.4170000000e-03, 1.3360000000e-02,
                5.0113000000e-02, 1.4621500000e-01, 3.0740500000e-01,
                2.6705100000e-01, -2.6213100000e-01, -6.3199900000e-01,
                -2.5160000000e-03, 6.0170100000e-01, 1.8204300000e-01},
      doubles_t{4.8911500000e+02, 1.4720800000e+02, 5.7127800000e+01,
                2.5042900000e+01, 1.1667200000e+01, 5.6825000000e+00,
                2.8132500000e+00, 1.3773300000e+00, 6.5653400000e-01,
                3.0070000000e-01, 1.3018700000e-01, 5.1582000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.4300000000e-04, -3.3700000000e-03, -1.8245000000e-02,
                -7.0772000000e-02, -2.1786200000e-01, -4.7743300000e-01,
                1.1012500000e-01, 1.0166160000e+00, -3.8078300000e-01,
                -8.4750000000e-01, 5.9344100000e-01, 2.7945300000e-01},
      doubles_t{4.8911500000e+02, 1.4720800000e+02, 5.7127800000e+01,
                2.5042900000e+01, 1.1667200000e+01, 5.6825000000e+00,
                2.8132500000e+00, 1.3773300000e+00, 6.5653400000e-01,
                3.0070000000e-01, 1.3018700000e-01, 5.1582000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-4.8000000000e-04, -4.1140000000e-03, -2.5128000000e-02,
                -9.0876000000e-02, -3.5693000000e-01, -4.9830400000e-01,
                1.2870220000e+00, -1.0200500000e-01, -1.4530250000e+00,
                1.4203760000e+00, -2.1296400000e-01, -5.0082400000e-01},
      doubles_t{4.8911500000e+02, 1.4720800000e+02, 5.7127800000e+01,
                2.5042900000e+01, 1.1667200000e+01, 5.6825000000e+00,
                2.8132500000e+00, 1.3773300000e+00, 6.5653400000e-01,
                3.0070000000e-01, 1.3018700000e-01, 5.1582000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1582000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0440000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.5384000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6229000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1260000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1750000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0666000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8324000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9704000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6570000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4539000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2600000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3156000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.2862000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1502000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4049900000e+00}));
    return abs_t(name, 26, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv5z_26

} // namespace chemcache
