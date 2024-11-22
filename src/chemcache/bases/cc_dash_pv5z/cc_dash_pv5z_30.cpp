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

abs_t cc_dash_pv5z_30() {
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
      doubles_t{1.0090140000e-07,  6.7047470000e-07,  3.8631660000e-06,
                1.7561080000e-05,  6.7471460000e-05,  2.2656520000e-04,
                6.8645960000e-04,  1.9195560000e-03,  5.0350360000e-03,
                1.2491830000e-02,  2.9311740000e-02,  6.4273670000e-02,
                1.2802920000e-01,  2.1961430000e-01,  2.9432850000e-01,
                2.5926390000e-01,  1.1541330000e-01,  1.7481790000e-02,
                1.2344700000e-03,  1.0308870000e-05,  5.4073860000e-05,
                -9.0209250000e-05, 1.2292050000e-05,  -1.2912990000e-05,
                5.0530070000e-06,  -4.0832950000e-06, 1.7740700000e-06,
                -4.5129620000e-07},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.1157410000e-08, -2.0703120000e-07, -1.1929350000e-06,
                -5.4227760000e-06, -2.0838260000e-05, -6.9989880000e-05,
                -2.1223190000e-04, -5.9441220000e-04, -1.5651210000e-03,
                -3.9128250000e-03, -9.3291740000e-03, -2.1099180000e-02,
                -4.4608570000e-02, -8.5438540000e-02, -1.3985180000e-01,
                -1.6749070000e-01, -7.8648970000e-02, 1.8714670000e-01,
                4.7737400000e-01,  3.9746080000e-01,  1.1191240000e-01,
                1.0210180000e-02,  1.5354660000e-03,  3.5132730000e-04,
                8.7003260000e-06,  8.0859280000e-06,  -4.4657790000e-06,
                1.2620790000e-06},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.1764680000e-08,  7.8186300000e-08,  4.5041710000e-07,
                2.0480430000e-06,  7.8675650000e-06,  2.6436790000e-05,
                8.0131270000e-05,  2.2461430000e-04,  5.9125140000e-04,
                1.4808200000e-03,  3.5339560000e-03,  8.0322400000e-03,
                1.7084550000e-02,  3.3235580000e-02,  5.5731780000e-02,
                7.0045260000e-02,  3.4608940000e-02,  -9.4174570000e-02,
                -3.1857850000e-01, -3.7401360000e-01, -7.2234240000e-03,
                4.6990510000e-01,  5.6500180000e-01,  2.1284460000e-01,
                1.0207730000e-02,  -3.2166470000e-03, 9.8139010000e-04,
                -2.1259530000e-04},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.2060660000e-08,  8.1101690000e-08,  4.6020500000e-07,
                2.1309590000e-06,  8.0045270000e-06,  2.7653820000e-05,
                8.0993100000e-05,  2.3672700000e-04,  5.9286090000e-04,
                1.5749590000e-03,  3.5171110000e-03,  8.6516050000e-03,
                1.6896280000e-02,  3.6868330000e-02,  5.4450030000e-02,
                8.6906710000e-02,  1.7404430000e-02,  -5.2125730000e-02,
                -5.5117260000e-01, -1.1543500000e+00, 2.1468370000e+00,
                2.4632420000e+00,  -6.0927930000e+00, 3.2604830000e+00,
                2.0784120000e+00,  -4.0424470000e+00, 1.7426680000e+00,
                1.5690180000e-01},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.7088660000e-09,  4.4166200000e-08,  2.5753650000e-07,
                1.1540120000e-06,  4.5135710000e-06,  1.4831940000e-05,
                4.6210340000e-05,  1.2526040000e-04,  3.4319940000e-04,
                8.2036100000e-04,  2.0668980000e-03,  4.4197890000e-03,
                1.0113670000e-02,  1.8094260000e-02,  3.4166620000e-02,
                3.5932650000e-02,  3.0741120000e-02,  -7.7526290000e-02,
                -1.6791940000e-01, -3.5624990000e-01, 1.7769820000e-01,
                4.7491290000e-01,  8.5219090000e-01,  -7.7248730000e-01,
                -1.4319140000e+00, 4.9919460000e-01,  5.4995850000e-01,
                2.4138940000e-01},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.3564930000e-09, -1.5659390000e-08, -9.0222050000e-08,
                -4.1017830000e-07, -1.5759920000e-06, -5.2945140000e-06,
                -1.6052610000e-05, -4.4982590000e-05, -1.1846240000e-04,
                -2.9658950000e-04, -7.0836590000e-04, -1.6097910000e-03,
                -3.4296790000e-03, -6.6778300000e-03, -1.1248160000e-02,
                -1.4190670000e-02, -7.1197830000e-03, 1.9642330000e-02,
                6.8660330000e-02,  8.4516910000e-02,  -2.5585930000e-03,
                -1.2468730000e-01, -1.9817320000e-01, -1.6269830000e-01,
                1.2126770000e-01,  4.1476240000e-01,  4.8408720000e-01,
                1.4192220000e-01},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.3178840000e-09, -6.1556070000e-08, -3.5734160000e-07,
                -1.6098810000e-06, -6.2550930000e-06, -2.0724440000e-05,
                -6.3920440000e-05, -1.7542690000e-04, -4.7366610000e-04,
                -1.1520850000e-03, -2.8465040000e-03, -6.2300270000e-03,
                -1.3898480000e-02, -2.5724960000e-02, -4.6749720000e-02,
                -5.3041360000e-02, -3.8466760000e-02, 9.7003090000e-02,
                2.7501320000e-01,  6.1491350000e-01,  -4.6482420000e-01,
                -1.4688170000e+00, -6.1233860000e-02, 2.5333260000e+00,
                -1.4838310000e+00, -9.4402430000e-01, 9.9343430000e-01,
                2.8341310000e-01},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.2880130000e-08,  7.9634940000e-08,  5.0283110000e-07,
                2.0448700000e-06,  8.9974620000e-06,  2.5480890000e-05,
                9.5025630000e-05,  2.0575900000e-04,  7.3222200000e-04,
                1.2790620000e-03,  4.5837980000e-03,  6.4884350000e-03,
                2.3696040000e-02,  2.3608410000e-02,  9.1869540000e-02,
                1.5845930000e-02,  1.8322090000e-01,  -3.6208580000e-01,
                -2.5620350000e-01, -2.6784160000e+00, 8.3186980000e+00,
                -7.7549970000e+00, 5.7989220000e-01,  4.4806790000e+00,
                -7.9856680000e+00, 7.9331330000e+00,  -2.6833390000e+00,
                -4.2356200000e-02},
      doubles_t{1.8741300000e+08, 3.3512700000e+07, 7.5047500000e+06,
                2.0178700000e+06, 6.2857800000e+05, 2.2016300000e+05,
                8.4595200000e+04, 3.4955500000e+04, 1.5290300000e+04,
                6.9949300000e+03, 3.3165900000e+03, 1.6192400000e+03,
                8.1034300000e+02, 4.1439300000e+02, 2.1605800000e+02,
                1.1463500000e+02, 6.1762300000e+01, 3.3692400000e+01,
                1.8531500000e+01, 1.0214400000e+01, 5.5942800000e+00,
                3.0096800000e+00, 1.5668400000e+00, 7.8143700000e-01,
                2.5811800000e-01, 1.3855600000e-01, 6.6171000000e-02,
                3.0986000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0986000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06,  2.5000000000e-05, 1.5400000000e-04,
                7.1100000000e-04,  2.7260000000e-03, 9.0780000000e-03,
                2.6647000000e-02,  6.8349000000e-02, 1.4622900000e-01,
                2.4636600000e-01,  3.0935700000e-01, 2.5183400000e-01,
                1.0757700000e-01,  1.8623000000e-02, 1.8400000000e-03,
                5.2600000000e-04,  1.3400000000e-04, 1.8000000000e-05,
                -3.0000000000e-06, 1.0000000000e-06},
      doubles_t{1.1034900000e+05, 2.6971300000e+04, 8.5600700000e+03,
                3.2488000000e+03, 1.3894100000e+03, 6.4304300000e+02,
                3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
                4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01,
                7.8326700000e+00, 4.2841000000e+00, 2.2890100000e+00,
                1.1943600000e+00, 5.9425800000e-01, 2.2691400000e-01,
                9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0000000000e-06, -9.0000000000e-06, -5.7000000000e-05,
                -2.6600000000e-04, -1.0230000000e-03, -3.4260000000e-03,
                -1.0178000000e-02, -2.6684000000e-02, -5.9044000000e-02,
                -1.0364500000e-01, -1.3718600000e-01, -9.8030000000e-02,
                7.8284000000e-02,  3.0689600000e-01,  4.0209100000e-01,
                2.8399600000e-01,  8.6928000000e-02,  5.8540000000e-03,
                -3.7000000000e-05, 1.0500000000e-04},
      doubles_t{1.1034900000e+05, 2.6971300000e+04, 8.5600700000e+03,
                3.2488000000e+03, 1.3894100000e+03, 6.4304300000e+02,
                3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
                4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01,
                7.8326700000e+00, 4.2841000000e+00, 2.2890100000e+00,
                1.1943600000e+00, 5.9425800000e-01, 2.2691400000e-01,
                9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0000000000e-06, -7.0000000000e-06, -4.4000000000e-05,
                -1.8700000000e-04, -7.8600000000e-04, -2.4190000000e-03,
                -7.8320000000e-03, -1.8892000000e-02, -4.6095000000e-02,
                -7.2187000000e-02, -1.1870600000e-01, -6.8845000000e-02,
                6.1255000000e-02,  5.4829400000e-01,  2.9728200000e-01,
                -2.7160100000e-01, -1.0379400000e+00, 2.8832900000e-01,
                5.8292400000e-01,  -1.0259000000e-02},
      doubles_t{1.1034900000e+05, 2.6971300000e+04, 8.5600700000e+03,
                3.2488000000e+03, 1.3894100000e+03, 6.4304300000e+02,
                3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
                4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01,
                7.8326700000e+00, 4.2841000000e+00, 2.2890100000e+00,
                1.1943600000e+00, 5.9425800000e-01, 2.2691400000e-01,
                9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        3.0000000000e-06, 1.9000000000e-05, 8.4000000000e-05, 3.3500000000e-04,
        1.0780000000e-03, 3.3420000000e-03, 8.4110000000e-03, 1.9569000000e-02,
        3.2567000000e-02, 4.7782000000e-02, 2.7957000000e-02, -2.4530000000e-02,
        -1.4649000000e-01, -1.4925700000e-01, -1.3631800000e-01,
        1.9895200000e-01, 6.5205700000e-01, 3.1011900000e-01, 5.5280000000e-03},
      doubles_t{
        2.6971300000e+04, 8.5600700000e+03, 3.2488000000e+03, 1.3894100000e+03,
        6.4304300000e+02, 3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
        4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01, 7.8326700000e+00,
        4.2841000000e+00, 2.2890100000e+00, 1.1943600000e+00, 5.9425800000e-01,
        2.2691400000e-01, 9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  1.1000000000e-05,  7.5000000000e-05,
                3.0700000000e-04,  1.3430000000e-03,  3.9790000000e-03,
                1.3508000000e-02,  3.1518000000e-02,  8.1224000000e-02,
                1.1774400000e-01,  2.5453700000e-01,  2.8046600000e-01,
                -1.3817560000e+00, -5.8450200000e-01, 2.8319550000e+00,
                -1.4030820000e+00, -1.0866400000e+00, 1.8537790000e+00,
                -1.1119850000e+00, 9.9530000000e-03},
      doubles_t{1.1034900000e+05, 2.6971300000e+04, 8.5600700000e+03,
                3.2488000000e+03, 1.3894100000e+03, 6.4304300000e+02,
                3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
                4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01,
                7.8326700000e+00, 4.2841000000e+00, 2.2890100000e+00,
                1.1943600000e+00, 5.9425800000e-01, 2.2691400000e-01,
                9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  9.0000000000e-06,  4.7000000000e-05,
                2.5800000000e-04,  8.4000000000e-04,  3.3290000000e-03,
                8.4050000000e-03,  2.6210000000e-02,  4.8326000000e-02,
                1.0642600000e-01,  1.1124300000e-01,  2.1963400000e-01,
                -3.7850600000e-01, -8.3428300000e-01, 3.4651000000e-02,
                1.9667590000e+00,  -1.2248450000e+00, -7.3773200000e-01,
                9.9164400000e-01,  -3.0089000000e-02},
      doubles_t{1.1034900000e+05, 2.6971300000e+04, 8.5600700000e+03,
                3.2488000000e+03, 1.3894100000e+03, 6.4304300000e+02,
                3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
                4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01,
                7.8326700000e+00, 4.2841000000e+00, 2.2890100000e+00,
                1.1943600000e+00, 5.9425800000e-01, 2.2691400000e-01,
                9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        2.0000000000e-06, 1.1000000000e-05, 4.9000000000e-05, 1.9000000000e-04,
        6.2800000000e-04, 1.8980000000e-03, 4.9050000000e-03, 1.1091000000e-02,
        1.9110000000e-02, 2.6404000000e-02, 1.7048000000e-02, -1.6777000000e-02,
        -7.2281000000e-02, -8.8834000000e-02, -8.5501000000e-02,
        5.1477000000e-02, 3.4826300000e-01, 5.2071300000e-01, 2.4111000000e-01},
      doubles_t{
        2.6971300000e+04, 8.5600700000e+03, 3.2488000000e+03, 1.3894100000e+03,
        6.4304300000e+02, 3.1400500000e+02, 1.5948900000e+02, 8.3666200000e+01,
        4.5157400000e+01, 2.4963000000e+01, 1.4004300000e+01, 7.8326700000e+00,
        4.2841000000e+00, 2.2890100000e+00, 1.1943600000e+00, 5.9425800000e-01,
        2.2691400000e-01, 9.0179000000e-02, 3.5718000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5718000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.1000000000e-04, 1.0600000000e-03, 6.1270000000e-03,
                2.4227000000e-02, 7.0568000000e-02, 1.5542200000e-01,
                2.4646500000e-01, 2.9085400000e-01, 2.6962300000e-01,
                1.9320200000e-01, 9.2824000000e-02, 1.9273000000e-02},
      doubles_t{8.5412230000e+02, 2.5794960000e+02, 1.0023190000e+02,
                4.4176210000e+01, 2.0784760000e+01, 1.0227560000e+01,
                5.1512410000e+00, 2.5908470000e+00, 1.2858010000e+00,
                6.2534400000e-01, 2.9394300000e-01, 1.2899700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.5400000000e-04, 1.4840000000e-03, 8.6150000000e-03,
                3.4328000000e-02, 1.0265600000e-01, 2.3462300000e-01,
                3.3085800000e-01, 1.6641400000e-01, -2.1915000000e-01,
                -4.1717900000e-01, -2.8798800000e-01, -7.6454000000e-02},
      doubles_t{8.5412230000e+02, 2.5794960000e+02, 1.0023190000e+02,
                4.4176210000e+01, 2.0784760000e+01, 1.0227560000e+01,
                5.1512410000e+00, 2.5908470000e+00, 1.2858010000e+00,
                6.2534400000e-01, 2.9394300000e-01, 1.2899700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.6800000000e-04, 1.6860000000e-03, 9.5520000000e-03,
                3.9169000000e-02, 1.1854000000e-01, 3.0601700000e-01,
                2.8003300000e-01, -3.5704400000e-01, -5.5669300000e-01,
                -3.4070000000e-03, 5.5504800000e-01, 3.0779800000e-01},
      doubles_t{8.5412230000e+02, 2.5794960000e+02, 1.0023190000e+02,
                4.4176210000e+01, 2.0784760000e+01, 1.0227560000e+01,
                5.1512410000e+00, 2.5908470000e+00, 1.2858010000e+00,
                6.2534400000e-01, 2.9394300000e-01, 1.2899700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.0200000000e-04, 2.2920000000e-03, 1.1948000000e-02,
                5.3489000000e-02, 1.5995700000e-01, 4.6332900000e-01,
                1.8920000000e-03, -9.5924800000e-01, 1.0855300000e-01,
                9.8188500000e-01, -3.0198700000e-01, -5.9589300000e-01},
      doubles_t{8.5412230000e+02, 2.5794960000e+02, 1.0023190000e+02,
                4.4176210000e+01, 2.0784760000e+01, 1.0227560000e+01,
                5.1512410000e+00, 2.5908470000e+00, 1.2858010000e+00,
                6.2534400000e-01, 2.9394300000e-01, 1.2899700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.2100000000e-04, -2.8130000000e-03, -1.8028000000e-02,
                -6.7985000000e-02, -2.8146800000e-01, -6.3964700000e-01,
                1.0829690000e+00, 4.4146900000e-01, -1.6892300000e+00,
                9.0111200000e-01, 6.1214400000e-01, -8.5295700000e-01},
      doubles_t{8.5412230000e+02, 2.5794960000e+02, 1.0023190000e+02,
                4.4176210000e+01, 2.0784760000e+01, 1.0227560000e+01,
                5.1512410000e+00, 2.5908470000e+00, 1.2858010000e+00,
                6.2534400000e-01, 2.9394300000e-01, 1.2899700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2899700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2166900000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3915000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5850000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7210000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{9.6299000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3924000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1951000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.4738000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4250000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5487000000e+00}));
    return abs_t(name, 30, r0, shells.begin(), shells.end());
} // cc_dash_pv5z_30

} // namespace chemcache
