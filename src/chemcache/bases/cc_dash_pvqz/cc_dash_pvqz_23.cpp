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

abs_t cc_dash_pvqz_23() {
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
      doubles_t{1.9662910000e-06, 1.5289670000e-05,  8.0439320000e-05,
                3.3976170000e-04, 1.2374050000e-03,  4.0257400000e-03,
                1.1919120000e-02, 3.2154310000e-02,  7.7854520000e-02,
                1.6163670000e-01, 2.6341940000e-01,  2.7999860000e-01,
                1.4734740000e-01, 6.8264050000e-02,  1.1944360000e-01,
                8.5075820000e-02, 1.1258140000e-02,  -3.3535210000e-04,
                1.0705710000e-04, -2.9796220000e-05, 4.6728160000e-05,
                -3.3212730000e-08},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1009230000e-06, -8.5632500000e-06, -4.5036780000e-05,
                -1.9038080000e-04, -6.9342840000e-04, -2.2624290000e-03,
                -6.7249910000e-03, -1.8373250000e-02, -4.5552600000e-02,
                -9.9969690000e-02, -1.8191110000e-01, -2.4587440000e-01,
                -1.6380890000e-01, 1.7958430000e-01,  5.3895470000e-01,
                3.7490980000e-01,  5.2002690000e-02,  -3.7758730000e-03,
                1.6868880000e-04,  -3.0638020000e-04, 3.2380960000e-04,
                -3.6722440000e-05},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.2788540000e-07,  1.7739420000e-06,  9.3199920000e-06,
                3.9464490000e-05,  1.4357200000e-04,  4.7029110000e-04,
                1.3999980000e-03,  3.8741110000e-03,  9.7659420000e-03,
                2.2441790000e-02,  4.3900910000e-02,  6.8434620000e-02,
                5.0965000000e-02,  -8.5359400000e-02, -3.4748040000e-01,
                -3.6723680000e-01, 2.3854080000e-01,  6.8881780000e-01,
                3.5099190000e-01,  1.8261440000e-02,  -5.5411340000e-03,
                2.7130020000e-03},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.2759480000e-08, -4.1071290000e-07, -2.1577110000e-06,
                -9.1372200000e-06, -3.3238850000e-05, -1.0889250000e-04,
                -3.2414070000e-04, -8.9723190000e-04, -2.2621050000e-03,
                -5.2030680000e-03, -1.0191230000e-02, -1.5948200000e-02,
                -1.1935880000e-02, 2.0183440000e-02,  8.6100710000e-02,
                9.4879170000e-02,  -7.3097530000e-02, -2.3496320000e-01,
                -2.5002970000e-01, 3.0260960000e-01,  6.0616350000e-01,
                2.5835740000e-01},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0034790000e-07, -7.9970610000e-07, -4.0603970000e-06,
                -1.7980700000e-05, -6.1813700000e-05, -2.1689150000e-04,
                -5.9424170000e-04, -1.8145060000e-03, -4.0648070000e-03,
                -1.0792230000e-02, -1.7522640000e-02, -3.5738730000e-02,
                -1.2630790000e-02, 1.3934940000e-02,  2.3050460000e-01,
                9.3657400000e-02,  7.8289600000e-03,  -1.0538540000e+00,
                3.8655910000e-01,  1.9257830000e+00,  -1.1464290000e+00,
                -5.7990930000e-01},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.5921480000e-07, -1.2783050000e-06, -6.4200970000e-06,
                -2.8836330000e-05, -9.7358900000e-05, -3.4914020000e-04,
                -9.3161940000e-04, -2.9349590000e-03, -6.3337030000e-03,
                -1.7606980000e-02, -2.6968480000e-02, -5.9820790000e-02,
                -1.5456130000e-02, 1.1030010000e-02,  4.2340060000e-01,
                1.4540880000e-01,  -2.1971400000e-01, -2.1943290000e+00,
                2.8574410000e+00,  -5.8324750000e-01, -1.6928450000e+00,
                1.6427260000e+00},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.8851990000e-07, -1.3945500000e-06, -7.8815190000e-06,
                -3.0279810000e-05, -1.2431600000e-04, -3.5072310000e-04,
                -1.2464400000e-03, -2.7874200000e-03, -9.0456180000e-03,
                -1.5229670000e-02, -4.4374960000e-02, -3.7935070000e-02,
                -8.5642160000e-02, 1.7156350000e-01,  2.1923740000e-01,
                8.1903930000e-01,  -2.4396980000e+00, 9.0486850000e-01,
                1.7069230000e+00,  -4.6758140000e+00, 6.0487340000e+00,
                -2.8156970000e+00},
      doubles_t{
        1.0251780000e+07, 1.5349200000e+06, 3.4930090000e+05, 9.8942050000e+04,
        3.2281360000e+04, 1.1655270000e+04, 4.5465480000e+03, 1.8861850000e+03,
        8.2294170000e+02, 3.7434390000e+02, 1.7630690000e+02, 8.5406610000e+01,
        4.1904350000e+01, 1.9155250000e+01, 9.5716170000e+00, 4.7800860000e+00,
        2.1618600000e+00, 1.0427110000e+00, 4.8028200000e-01, 1.1763800000e-01,
        5.6174000000e-02, 2.5661000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.1000000000e-05, 1.0100000000e-04, 5.9200000000e-04,
                2.6900000000e-03, 1.0063000000e-02, 3.1654000000e-02,
                8.3478000000e-02, 1.7722600000e-01, 2.8839800000e-01,
                3.2943200000e-01, 2.1213000000e-01, 5.5882000000e-02,
                2.0970000000e-03, -1.4710000000e-03, -9.4800000000e-04,
                -5.5000000000e-05, -1.3000000000e-05, 2.0000000000e-06},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-4.0000000000e-06, -3.4000000000e-05, -2.0200000000e-04,
                -9.2000000000e-04, -3.4590000000e-03, -1.1025000000e-02,
                -2.9732000000e-02, -6.5349000000e-02, -1.1052700000e-01,
                -1.3453700000e-01, -5.2437000000e-02, 1.8221600000e-01,
                3.9801800000e-01, 3.8915300000e-01, 1.7271700000e-01,
                1.8746000000e-02, -5.9400000000e-04, 4.2200000000e-04},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 1.0000000000e-05, 6.0000000000e-05,
                2.7500000000e-04, 1.0330000000e-03, 3.2960000000e-03,
                8.8960000000e-03, 1.9613000000e-02, 3.3318000000e-02,
                4.1009000000e-02, 1.4959000000e-02, -6.4256000000e-02,
                -1.4181700000e-01, -1.7396600000e-01, -4.1132000000e-02,
                5.6399500000e-01, 5.2389600000e-01, 3.5689000000e-02},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06, 2.0000000000e-05, 1.2900000000e-04,
                5.3600000000e-04, 2.2030000000e-03, 6.4280000000e-03,
                1.9071000000e-02, 3.7712000000e-02, 7.4267000000e-02,
                7.5260000000e-02, 6.0998000000e-02, -2.0097600000e-01,
                -3.4534300000e-01, -5.0183500000e-01, 1.1436310000e+00,
                5.0448700000e-01, -9.4410300000e-01, -7.5778000000e-02},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.0000000000e-06, -2.9000000000e-05, -1.4100000000e-04,
                -7.7300000000e-04, -2.4220000000e-03, -9.2670000000e-03,
                -2.0555000000e-02, -5.6485000000e-02, -7.1493000000e-02,
                -1.5779100000e-01, -2.1552000000e-02, 2.3795000000e-01,
                1.1318370000e+00, -1.0891800000e+00, -1.1257730000e+00,
                2.1681890000e+00, -1.1102630000e+00, -1.8401000000e-01},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 8.0000000000e-06, 4.9000000000e-05,
                2.2100000000e-04, 8.3800000000e-04, 2.6560000000e-03,
                7.2330000000e-03, 1.5808000000e-02, 2.7172000000e-02,
                3.2801000000e-02, 1.2338000000e-02, -5.5418000000e-02,
                -1.1560400000e-01, -1.3909200000e-01, 6.8200000000e-04,
                3.5200900000e-01, 5.4220500000e-01, 2.3135800000e-01},
      doubles_t{2.8563320000e+04, 6.7571910000e+03, 2.1946520000e+03,
                8.4052570000e+02, 3.5766660000e+02, 1.6375340000e+02,
                7.9144660000e+01, 3.9838580000e+01, 2.0585360000e+01,
                1.0843070000e+01, 5.7779690000e+00, 3.0431700000e+00,
                1.5702390000e+00, 7.9497500000e-01, 3.8898000000e-01,
                1.5732700000e-01, 6.6566000000e-02, 2.7861000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7861000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.5000000000e-04, 3.1300000000e-03, 1.5478000000e-02,
                           5.0829000000e-02, 1.2754800000e-01, 2.3044700000e-01,
                           3.0405800000e-01, 3.0836600000e-01, 2.2902200000e-01,
                           1.0136000000e-01, 1.2961000000e-02},
                 doubles_t{2.3800900000e+02, 7.1189600000e+01, 2.7394600000e+01,
                           1.1717300000e+01, 5.3248000000e+00, 2.5169800000e+00,
                           1.1907500000e+00, 5.5180900000e-01, 2.4692600000e-01,
                           1.0510000000e-01, 4.1393000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.7000000000e-04, -3.3170000000e-03, -1.6474000000e-02,
                -5.4839000000e-02, -1.4083900000e-01, -2.3955300000e-01,
                -2.5191000000e-01, -2.6106000000e-02, 3.8777000000e-01,
                4.8672400000e-01, 1.3422400000e-01},
      doubles_t{2.3800900000e+02, 7.1189600000e+01, 2.7394600000e+01,
                1.1717300000e+01, 5.3248000000e+00, 2.5169800000e+00,
                1.1907500000e+00, 5.5180900000e-01, 2.4692600000e-01,
                1.0510000000e-01, 4.1393000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-5.7500000000e-04, -5.1570000000e-03, -2.5929000000e-02,
                -8.8059000000e-02, -2.3117900000e-01, -3.3835200000e-01,
                -8.5030000000e-02, 6.2679900000e-01, 3.8347300000e-01,
                -6.1026200000e-01, -2.8410400000e-01},
      doubles_t{2.3800900000e+02, 7.1189600000e+01, 2.7394600000e+01,
                1.1717300000e+01, 5.3248000000e+00, 2.5169800000e+00,
                1.1907500000e+00, 5.5180900000e-01, 2.4692600000e-01,
                1.0510000000e-01, 4.1393000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{7.3400000000e-04, 7.2580000000e-03, 3.3958000000e-02,
                1.2985500000e-01, 3.3164800000e-01, 3.5698300000e-01,
                -6.9343900000e-01, -5.7685000000e-01, 1.2280290000e+00,
                -3.1674100000e-01, -4.9428400000e-01},
      doubles_t{2.3800900000e+02, 7.1189600000e+01, 2.7394600000e+01,
                1.1717300000e+01, 5.3248000000e+00, 2.5169800000e+00,
                1.1907500000e+00, 5.5180900000e-01, 2.4692600000e-01,
                1.0510000000e-01, 4.1393000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1393000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0648000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.8400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1538000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.6450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3165000000e+00}));
    return abs_t(name, 23, r0, shells.begin(), shells.end());
} // cc_dash_pvqz_23

} // namespace chemcache
