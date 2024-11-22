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

abs_t aug_dash_cc_dash_pv5z_27() {
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
      doubles_t{8.6621830000e-08,  5.7717630000e-07,  3.3387770000e-06,
                1.5256390000e-05,  5.8964290000e-05,  1.9927020000e-04,
                6.0777890000e-04,  1.7109610000e-03,  4.5174740000e-03,
                1.1277820000e-02,  2.6612270000e-02,  5.8606370000e-02,
                1.1694820000e-01,  1.9997460000e-01,  2.6437840000e-01,
                2.2493410000e-01,  9.9062400000e-02,  5.0714820000e-02,
                8.9026850000e-02,  7.5861520000e-02,  2.2142220000e-02,
                1.7175310000e-03,  -5.4094070000e-06, -6.7184330000e-05,
                -7.1786000000e-07, 6.0325600000e-06,  6.1263130000e-06,
                2.5091120000e-06},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.5561830000e-08, -3.0358470000e-07, -1.7561660000e-06,
                -8.0248210000e-06, -3.1017680000e-05, -1.0484280000e-04,
                -3.1993140000e-04, -9.0161050000e-04, -2.3863900000e-03,
                -5.9881000000e-03, -1.4280630000e-02, -3.2117320000e-02,
                -6.6784960000e-02, -1.2366230000e-01, -1.9080180000e-01,
                -2.1154680000e-01, -9.9222280000e-02, 1.7229460000e-01,
                4.6192150000e-01,  3.9962470000e-01,  1.1627980000e-01,
                8.6044980000e-03,  -1.2652930000e-03, -7.6485860000e-04,
                -6.4043520000e-05, 6.3726320000e-05,  1.6629650000e-05,
                1.5700130000e-05},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0415620000e-08,  6.9414610000e-08,  4.0144700000e-07,
                1.8350090000e-06,  7.0904290000e-06,  2.3980840000e-05,
                7.3162860000e-05,  2.0649440000e-04,  5.4709300000e-04,
                1.3793590000e-03,  3.3115030000e-03,  7.5707740000e-03,
                1.6173810000e-02,  3.1582990000e-02,  5.3071110000e-02,
                6.6893190000e-02,  3.3388520000e-02,  -8.5532130000e-02,
                -2.9832110000e-01, -3.7206810000e-01, -3.3132810000e-02,
                4.4767130000e-01,  5.7112230000e-01,  2.3463700000e-01,
                1.5058340000e-02,  -3.4443530000e-03, 3.2123690000e-03,
                1.6750180000e-04},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.2007970000e-09, -1.4668500000e-08, -8.4822590000e-08,
                -3.8777840000e-07, -1.4981050000e-06, -5.0679130000e-06,
                -1.5457630000e-05, -4.3642200000e-05, -1.1558670000e-04,
                -2.9157330000e-04, -6.9973880000e-04, -1.6011430000e-03,
                -3.4202340000e-03, -6.6918470000e-03, -1.1258250000e-02,
                -1.4275620000e-02, -7.1347680000e-03, 1.8483410000e-02,
                6.7244840000e-02,  8.7043560000e-02,  4.9980140000e-03,
                -1.2385560000e-01, -1.9895200000e-01, -1.7599400000e-01,
                7.9445200000e-02,  4.0249330000e-01,  5.0340630000e-01,
                1.7304810000e-01},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.6002400000e-09, -5.0060640000e-08, -2.9389040000e-07,
                -1.3193570000e-06, -5.2117100000e-06, -1.7152610000e-05,
                -5.4109680000e-05, -1.4664550000e-04, -4.0772400000e-04,
                -9.7191380000e-04, -2.4896390000e-03, -5.2910210000e-03,
                -1.2330310000e-02, -2.1779190000e-02, -4.2153060000e-02,
                -4.2817630000e-02, -4.0648180000e-02, 9.6662740000e-02,
                1.9524510000e-01,  5.1394520000e-01,  -2.5450050000e-01,
                -9.2195390000e-01, -1.0086430000e+00, 2.2646620000e+00,
                5.8034640000e-01,  -2.1053590000e+00, 8.9094440000e-02,
                9.6257620000e-01},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.3057520000e-09, -2.8483850000e-08, -1.6629770000e-07,
                -7.5154590000e-07, -2.9446840000e-06, -9.7895640000e-06,
                -3.0503780000e-05, -8.3919090000e-05, -2.2920770000e-04,
                -5.5782240000e-04, -1.3950520000e-03, -3.0459220000e-03,
                -6.8734440000e-03, -1.2597700000e-02, -2.3148660000e-02,
                -2.5483970000e-02, -1.9580460000e-02, 4.7728970000e-02,
                1.1492970000e-01,  2.2570420000e-01,  -6.5541480000e-02,
                -2.4934250000e-01, -6.6863890000e-01, 2.0411080000e-01,
                1.2279070000e+00,  1.8575180000e-01,  -4.7632810000e-01,
                -6.7048140000e-01},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0705340000e-08, -7.0387140000e-08, -4.1419230000e-07,
                -1.8540490000e-06, -7.3500410000e-06, -2.4084110000e-05,
                -7.6387310000e-05, -2.0569450000e-04, -5.7637210000e-04,
                -1.3625150000e-03, -3.5274240000e-03, -7.4242930000e-03,
                -1.7582550000e-02, -3.0750950000e-02, -6.1332170000e-02,
                -6.1295590000e-02, -6.5436930000e-02, 1.5239380000e-01,
                3.4849880000e-01,  1.6044000000e+00,  -3.2744650000e+00,
                -1.4492520000e+00, 7.5471700000e+00,  -6.9914700000e+00,
                3.1552020000e+00,  9.7841410000e-01,  -3.6433640000e+00,
                1.9481030000e+00},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.1970160000e-09,  6.1870330000e-08,  3.5353270000e-07,
                1.6395540000e-06,  6.2235940000e-06,  2.1514910000e-05,
                6.3897190000e-05,  1.8632100000e-04,  4.7496940000e-04,
                1.2529740000e-03,  2.8592090000e-03,  6.9406240000e-03,
                1.3895730000e-02,  2.9599720000e-02,  4.5163140000e-02,
                6.8418080000e-02,  1.8397130000e-02,  -5.1862290000e-02,
                -3.9745220000e-01, -6.6067750000e-01, 5.2155810000e-01,
                2.8381460000e+00,  -2.9706580000e+00, -8.5075570000e-01,
                4.2455590000e+00,  -2.9372930000e+00, -1.1409540000e+00,
                1.4381650000e+00},
      doubles_t{1.6103000000e+08, 2.8755800000e+07, 6.4232400000e+06,
                1.7212000000e+06, 5.3399400000e+05, 1.8619400000e+05,
                7.1201700000e+04, 2.9276200000e+04, 1.2742300000e+04,
                5.8005100000e+03, 2.7369600000e+03, 1.3299900000e+03,
                6.6261100000e+02, 3.3741500000e+02, 1.7523400000e+02,
                9.2643400000e+01, 4.9757600000e+01, 2.7072800000e+01,
                1.4861000000e+01, 8.1812100000e+00, 4.4794100000e+00,
                2.4119100000e+00, 1.2584400000e+00, 6.3123400000e-01,
                2.1529700000e-01, 1.1851000000e-01, 5.7209000000e-02,
                2.7166000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7166000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2900000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06, 2.1000000000e-05, 1.3000000000e-04,
                6.0500000000e-04, 2.3310000000e-03, 7.8070000000e-03,
                2.3089000000e-02, 5.9908000000e-02, 1.3112000000e-01,
                2.2832300000e-01, 3.0338900000e-01, 2.7353200000e-01,
                1.3399100000e-01, 2.5516000000e-02, -6.9000000000e-04,
                -1.8610000000e-03, -7.8800000000e-04, -4.8000000000e-05,
                -2.0000000000e-06},
      doubles_t{
        9.7512600000e+04, 2.3815100000e+04, 7.5404300000e+03, 2.8540800000e+03,
        1.2176900000e+03, 5.6251500000e+02, 2.7427400000e+02, 1.3911400000e+02,
        7.2859800000e+01, 3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01,
        6.7691200000e+00, 3.7052800000e+00, 1.9790100000e+00, 1.0312700000e+00,
        5.1360900000e-01, 2.0426400000e-01, 8.3878000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0000000000e-06, -7.0000000000e-06, -4.6000000000e-05,
                -2.1500000000e-04, -8.3200000000e-04, -2.8010000000e-03,
                -8.3730000000e-03, -2.2160000000e-02, -5.0043000000e-02,
                -9.0502000000e-02, -1.2666800000e-01, -1.1131900000e-01,
                3.2705000000e-02,  2.6418900000e-01,  4.0160300000e-01,
                3.2442200000e-01,  1.2131800000e-01,  1.0815000000e-02,
                -2.5100000000e-04, 2.2700000000e-04},
      doubles_t{9.7512600000e+04, 2.3815100000e+04, 7.5404300000e+03,
                2.8540800000e+03, 1.2176900000e+03, 5.6251500000e+02,
                2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
                3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01,
                6.7691200000e+00, 3.7052800000e+00, 1.9790100000e+00,
                1.0312700000e+00, 5.1360900000e-01, 2.0426400000e-01,
                8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  1.1000000000e-05,  5.4000000000e-05,
                3.0700000000e-04,  9.7200000000e-04,  4.0050000000e-03,
                9.9110000000e-03,  3.2303000000e-02,  5.9473000000e-02,
                1.4088900000e-01,  1.3932200000e-01,  4.5638400000e-01,
                -9.8631300000e-01, -1.4009310000e+00, 2.8281690000e+00,
                -4.5441900000e-01, -1.9262680000e+00, 2.0851170000e+00,
                -9.6537700000e-01, -1.0770400000e-01},
      doubles_t{9.7512600000e+04, 2.3815100000e+04, 7.5404300000e+03,
                2.8540800000e+03, 1.2176900000e+03, 5.6251500000e+02,
                2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
                3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01,
                6.7691200000e+00, 3.7052800000e+00, 1.9790100000e+00,
                1.0312700000e+00, 5.1360900000e-01, 2.0426400000e-01,
                8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        3.0000000000e-06, 1.6000000000e-05, 7.3000000000e-05, 2.9000000000e-04,
        9.5100000000e-04, 2.9220000000e-03, 7.5360000000e-03, 1.7566000000e-02,
        3.0770000000e-02, 4.5854000000e-02, 3.6789000000e-02, -1.0786000000e-02,
        -1.2314900000e-01, -1.7051700000e-01, -1.5956300000e-01,
        1.3565100000e-01, 6.2549300000e-01, 3.7706700000e-01, 1.6199000000e-02},
      doubles_t{
        2.3815100000e+04, 7.5404300000e+03, 2.8540800000e+03, 1.2176900000e+03,
        5.6251500000e+02, 2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
        3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01, 6.7691200000e+00,
        3.7052800000e+00, 1.9790100000e+00, 1.0312700000e+00, 5.1360900000e-01,
        2.0426400000e-01, 8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  5.0000000000e-06,  3.5000000000e-05,
                1.4800000000e-04,  6.3400000000e-04,  1.9350000000e-03,
                6.3830000000e-03,  1.5340000000e-02,  3.8654000000e-02,
                6.1506000000e-02,  1.0635400000e-01,  7.2659000000e-02,
                2.6700000000e-04,  -4.2278600000e-01, -4.2361700000e-01,
                6.1196000000e-02,  1.1410350000e+00,  -7.0768000000e-02,
                -7.1971100000e-01, -1.8515000000e-02},
      doubles_t{9.7512600000e+04, 2.3815100000e+04, 7.5404300000e+03,
                2.8540800000e+03, 1.2176900000e+03, 5.6251500000e+02,
                2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
                3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01,
                6.7691200000e+00, 3.7052800000e+00, 1.9790100000e+00,
                1.0312700000e+00, 5.1360900000e-01, 2.0426400000e-01,
                8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06,  8.0000000000e-06,  3.8000000000e-05,
                2.1700000000e-04,  6.9000000000e-04,  2.8260000000e-03,
                6.9830000000e-03,  2.2541000000e-02,  4.1335000000e-02,
                9.6486000000e-02,  9.6311000000e-02,  2.1875800000e-01,
                -1.9916200000e-01, -7.5789600000e-01, -4.8819400000e-01,
                2.0656540000e+00,  -6.9018200000e-01, -1.2064820000e+00,
                1.0635780000e+00,  4.0615000000e-02},
      doubles_t{9.7512600000e+04, 2.3815100000e+04, 7.5404300000e+03,
                2.8540800000e+03, 1.2176900000e+03, 5.6251500000e+02,
                2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
                3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01,
                6.7691200000e+00, 3.7052800000e+00, 1.9790100000e+00,
                1.0312700000e+00, 5.1360900000e-01, 2.0426400000e-01,
                8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{
        2.0000000000e-06, 1.0000000000e-05, 4.7000000000e-05, 1.8500000000e-04,
        6.1300000000e-04, 1.8600000000e-03, 4.8640000000e-03, 1.1175000000e-02,
        1.9911000000e-02, 2.8809000000e-02, 2.4090000000e-02, -8.6610000000e-03,
        -7.2953000000e-02, -1.0667600000e-01, -1.0791200000e-01,
        3.9615000000e-02, 3.5604000000e-01, 5.2490900000e-01, 2.3089600000e-01},
      doubles_t{
        2.3815100000e+04, 7.5404300000e+03, 2.8540800000e+03, 1.2176900000e+03,
        5.6251500000e+02, 2.7427400000e+02, 1.3911400000e+02, 7.2859800000e+01,
        3.9244800000e+01, 2.1643000000e+01, 1.2114900000e+01, 6.7691200000e+00,
        3.7052800000e+00, 1.9790100000e+00, 1.0312700000e+00, 5.1360900000e-01,
        2.0426400000e-01, 8.3878000000e-02, 3.4142000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4142000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3900000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.9500000000e-04, 1.8200000000e-03, 9.8600000000e-03,
                3.5417000000e-02, 9.4882000000e-02, 1.8865100000e-01,
                2.7039500000e-01, 2.9634700000e-01, 2.5655900000e-01,
                1.6405500000e-01, 6.0567000000e-02, 6.0430000000e-03},
      doubles_t{4.9524800000e+02, 1.4903800000e+02, 5.7872800000e+01,
                2.5379600000e+01, 1.1859400000e+01, 5.7963000000e+00,
                2.8727500000e+00, 1.4078360000e+00, 6.7240700000e-01,
                3.0883500000e-01, 1.3429700000e-01, 5.3619000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-2.2600000000e-04, -2.1170000000e-03, -1.1506000000e-02,
                -4.1673000000e-02, -1.1485000000e-01, -2.3003600000e-01,
                -2.7872200000e-01, -1.0904400000e-01, 2.3091800000e-01,
                4.3417800000e-01, 3.1468100000e-01, 5.6026000000e-02},
      doubles_t{4.9524800000e+02, 1.4903800000e+02, 5.7872800000e+01,
                2.5379600000e+01, 1.1859400000e+01, 5.7963000000e+00,
                2.8727500000e+00, 1.4078360000e+00, 6.7240700000e-01,
                3.0883500000e-01, 1.3429700000e-01, 5.3619000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-4.1900000000e-04, -4.0790000000e-03, -2.1895000000e-02,
                -8.3031000000e-02, -2.5689100000e-01, -5.0375100000e-01,
                3.4976100000e-01, 9.3487000000e-01, -6.3839400000e-01,
                -6.3725600000e-01, 6.2545600000e-01, 2.1383500000e-01},
      doubles_t{4.9524800000e+02, 1.4903800000e+02, 5.7872800000e+01,
                2.5379600000e+01, 1.1859400000e+01, 5.7963000000e+00,
                2.8727500000e+00, 1.4078360000e+00, 6.7240700000e-01,
                3.0883500000e-01, 1.3429700000e-01, 5.3619000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.9800000000e-04, 2.8410000000e-03, 1.5380000000e-02,
                5.6932000000e-02, 1.6201800000e-01, 3.2512400000e-01,
                2.0714400000e-01, -3.6874000000e-01, -5.8084900000e-01,
                1.2670600000e-01, 5.8760600000e-01, 1.4627900000e-01},
      doubles_t{4.9524800000e+02, 1.4903800000e+02, 5.7872800000e+01,
                2.5379600000e+01, 1.1859400000e+01, 5.7963000000e+00,
                2.8727500000e+00, 1.4078360000e+00, 6.7240700000e-01,
                3.0883500000e-01, 1.3429700000e-01, 5.3619000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-6.5000000000e-04, -4.6330000000e-03, -3.2519000000e-02,
                -1.0107800000e-01, -4.6501900000e-01, -3.1860100000e-01,
                1.5149910000e+00, -7.6126700000e-01, -9.6716800000e-01,
                1.4433300000e+00, -4.7420000000e-01, -3.5861600000e-01},
      doubles_t{4.9524800000e+02, 1.4903800000e+02, 5.7872800000e+01,
                2.5379600000e+01, 1.1859400000e+01, 5.7963000000e+00,
                2.8727500000e+00, 1.4078360000e+00, 6.7240700000e-01,
                3.0883500000e-01, 1.3429700000e-01, 5.3619000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3619000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1410000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5920000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0187000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0606000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7260000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2655000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7723000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3190000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9748000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9511000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5513000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{6.3878000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6950000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7051400000e+00}));
    return abs_t(name, 27, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv5z_27

} // namespace chemcache
