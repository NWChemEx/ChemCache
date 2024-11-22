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

#include "aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_29() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(q+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9222100000e-06,  1.4946550000e-05,  7.8631060000e-05,
                3.3211890000e-04,  1.2096490000e-03,  3.9362310000e-03,
                1.1660310000e-02,  3.1495120000e-02,  7.6463850000e-02,
                1.5951620000e-01,  2.6189830000e-01,  2.8167080000e-01,
                1.5095870000e-01,  7.0419550000e-02,  1.1800390000e-01,
                8.0151870000e-02,  1.0201220000e-02,  -3.7204640000e-04,
                -2.2068860000e-05, -2.4010850000e-05, 1.5298160000e-05,
                -5.1511670000e-06},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.0826420000e-06, -8.4208690000e-06, -4.4286250000e-05,
                -1.8720710000e-04, -6.8192650000e-04, -2.2254850000e-03,
                -6.6196200000e-03, -1.8114000000e-02, -4.5066220000e-02,
                -9.9554430000e-02, -1.8305580000e-01, -2.5105770000e-01,
                -1.7034170000e-01, 1.9000770000e-01,  5.5070480000e-01,
                3.6425320000e-01,  4.8689770000e-02,  -3.7999930000e-03,
                -3.1563080000e-04, -2.4304110000e-04, 1.5536450000e-04,
                -5.1802560000e-05},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.3896050000e-07,  1.8600090000e-06,  9.7725800000e-06,
                4.1376190000e-05,  1.5055880000e-04,  4.9321750000e-04,
                1.4694670000e-03,  4.0714870000e-03,  1.0300910000e-02,
                2.3814080000e-02,  4.7070080000e-02,  7.4244500000e-02,
                5.6576570000e-02,  -9.7608620000e-02, -3.8563990000e-01,
                -3.5403120000e-01, 2.9839590000e-01,  6.8618090000e-01,
                3.1509890000e-01,  1.3945640000e-02,  -3.6116400000e-03,
                2.2891150000e-03},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.8234360000e-08, -3.7538580000e-07, -1.9727430000e-06,
                -8.3499610000e-06, -3.0395350000e-05, -9.9530100000e-05,
                -2.9671880000e-04, -8.2172930000e-04, -2.0813770000e-03,
                -4.8111490000e-03, -9.5376830000e-03, -1.5073510000e-02,
                -1.1610880000e-02, 2.0327160000e-02,  8.3274030000e-02,
                7.9531250000e-02,  -7.9015240000e-02, -2.0543120000e-01,
                -1.9851340000e-01, 2.6271330000e-01,  5.7940310000e-01,
                3.0541490000e-01},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.8493560000e-07, -1.4619050000e-06, -7.5105340000e-06,
                -3.2749890000e-05, -1.1482390000e-04, -3.9359030000e-04,
                -1.1106170000e-03, -3.2841740000e-03, -7.6972350000e-03,
                -1.9590130000e-02, -3.4447120000e-02, -6.5071760000e-02,
                -3.2694160000e-02, 4.7389150000e-02,  4.3337120000e-01,
                3.1895210000e-01,  -7.2819550000e-01, -1.6631160000e+00,
                2.5132670000e+00,  -3.7368610000e-01, -1.4791020000e+00,
                1.3158940000e+00},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-8.8766350000e-08, -7.0321350000e-07, -3.6013660000e-06,
                -1.5768710000e-05, -5.4997260000e-05, -1.8970550000e-04,
                -5.3116560000e-04, -1.5845930000e-03, -3.6713160000e-03,
                -9.4588930000e-03, -1.6301340000e-02, -3.1454020000e-02,
                -1.4447120000e-02, 2.0378560000e-02,  1.9973070000e-01,
                9.4759260000e-02,  -7.9353980000e-02, -7.5913530000e-01,
                2.3650790000e-01,  1.4226040000e+00,  -3.0662030000e-01,
                -9.6193100000e-01},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5303000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.9809510000e-07, -1.4461910000e-06, -8.3263080000e-06,
                -3.1193800000e-05, -1.3208420000e-04, -3.5853620000e-04,
                -1.3339960000e-03, -2.8247280000e-03, -9.8063330000e-03,
                -1.5293230000e-02, -4.9565730000e-02, -3.5834680000e-02,
                -1.0449060000e-01, 2.1058540000e-01,  1.6467320000e-01,
                1.1505130000e+00,  -3.4837320000e+00, 2.5166360000e+00,
                4.4290640000e-01,  -3.4043100000e+00, 4.6278500000e+00,
                -2.2401200000e+00},
      doubles_t{
        1.6665490000e+07, 2.4952130000e+06, 5.6785070000e+05, 1.6085310000e+05,
        5.2481860000e+04, 1.8948850000e+04, 7.3916550000e+03, 3.0665170000e+03,
        1.3379940000e+03, 6.0878300000e+02, 2.8689800000e+02, 1.3915300000e+02,
        6.8494550000e+01, 3.1601080000e+01, 1.5991580000e+01, 8.0565100000e+00,
        3.6903290000e+00, 1.7629020000e+00, 8.0126900000e-01, 1.8424300000e-01,
        8.1793000000e-02, 3.5303000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5240000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.1000000000e-05, 9.7000000000e-05, 5.6900000000e-04,
                2.5890000000e-03, 9.7420000000e-03, 3.0930000000e-02,
                8.2530000000e-02, 1.7775800000e-01, 2.9222500000e-01,
                3.3034900000e-01, 2.0686300000e-01, 5.2793000000e-02,
                1.7180000000e-03, -1.7010000000e-03, -9.7600000000e-04,
                -5.4000000000e-05, 1.0000000000e-07, -1.0000000000e-06},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-4.0000000000e-06, -3.5000000000e-05, -2.0600000000e-04,
                -9.4100000000e-04, -3.5560000000e-03, -1.1461000000e-02,
                -3.1309000000e-02, -7.0139000000e-02, -1.2056200000e-01,
                -1.4355200000e-01, -4.0327000000e-02, 2.0846600000e-01,
                4.0872400000e-01, 3.6820100000e-01, 1.5965300000e-01,
                1.5955000000e-02, -1.6810000000e-03, 6.3100000000e-04},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06, 1.3000000000e-05, 7.9000000000e-05,
                3.4900000000e-04, 1.3730000000e-03, 4.2610000000e-03,
                1.2146000000e-02, 2.6064000000e-02, 4.7907000000e-02,
                5.2829000000e-02, 2.1575000000e-02, -1.1778300000e-01,
                -1.8839700000e-01, -2.0191600000e-01, 1.7374300000e-01,
                6.3684800000e-01, 3.5151800000e-01, 1.2383000000e-02},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06, 2.5000000000e-05, 1.6300000000e-04,
                6.7600000000e-04, 2.8080000000e-03, 8.2550000000e-03,
                2.4946000000e-02, 5.0160000000e-02, 1.0028900000e-01,
                1.0231900000e-01, 7.7241000000e-02, -3.8444900000e-01,
                -4.5135400000e-01, -9.4778000000e-02, 1.1731180000e+00,
                -5.4133000000e-02, -7.2824400000e-01, 4.2700000000e-04},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0000000000e-06, 3.7000000000e-05, 1.9000000000e-04,
                9.7800000000e-04, 3.3080000000e-03, 1.1980000000e-02,
                2.9207000000e-02, 7.5122000000e-02, 1.0783500000e-01,
                2.0737900000e-01, 6.0406000000e-02, -8.9497700000e-01,
                -5.5840200000e-01, 1.8997060000e+00, -5.0222400000e-01,
                -1.1824130000e+00, 1.0445070000e+00, 1.2236000000e-02},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 7.0000000000e-06, 4.0000000000e-05,
                1.8200000000e-04, 6.9700000000e-04, 2.2200000000e-03,
                6.1630000000e-03, 1.3631000000e-02, 2.4033000000e-02,
                2.7862000000e-02, 7.9750000000e-03, -5.1133000000e-02,
                -9.2305000000e-02, -1.0645800000e-01, 3.7380000000e-03,
                3.1054100000e-01, 5.4190300000e-01, 2.8250200000e-01},
      doubles_t{4.8218410000e+04, 1.1420550000e+04, 3.7135690000e+03,
                1.4237620000e+03, 6.0660410000e+02, 2.7825640000e+02,
                1.3484490000e+02, 6.8122540000e+01, 3.5400880000e+01,
                1.8790290000e+01, 1.0083560000e+01, 5.3461200000e+00,
                2.7752280000e+00, 1.4075640000e+00, 6.8721300000e-01,
                2.4682200000e-01, 9.6220000000e-02, 3.7037000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7037000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4260000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.5600000000e-04, 2.3780000000e-03, 1.2722000000e-02,
                           4.5012000000e-02, 1.1719600000e-01, 2.1855500000e-01,
                           2.9144800000e-01, 2.9713800000e-01, 2.3609200000e-01,
                           1.3144400000e-01, 3.4345000000e-02},
                 doubles_t{5.2788600000e+02, 1.5873300000e+02, 6.1555600000e+01,
                           2.6883400000e+01, 1.2517500000e+01, 6.0697700000e+00,
                           2.9583190000e+00, 1.4141220000e+00, 6.5332900000e-01,
                           2.8669200000e-01, 1.1501400000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.1600000000e-04, -2.9580000000e-03, -1.5877000000e-02,
                -5.6864000000e-02, -1.5358400000e-01, -2.8113900000e-01,
                -2.6851800000e-01, 2.3130000000e-02, 3.5662000000e-01,
                4.1847000000e-01, 1.9282700000e-01},
      doubles_t{5.2788600000e+02, 1.5873300000e+02, 6.1555600000e+01,
                2.6883400000e+01, 1.2517500000e+01, 6.0697700000e+00,
                2.9583190000e+00, 1.4141220000e+00, 6.5332900000e-01,
                2.8669200000e-01, 1.1501400000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-5.9300000000e-04, -5.1060000000e-03, -2.9922000000e-02,
                -1.0329300000e-01, -3.5956800000e-01, -4.1280700000e-01,
                7.4768900000e-01, 5.4251500000e-01, -9.7770600000e-01,
                -8.8603000000e-02, 6.5507000000e-01},
      doubles_t{5.2788600000e+02, 1.5873300000e+02, 6.1555600000e+01,
                2.6883400000e+01, 1.2517500000e+01, 6.0697700000e+00,
                2.9583190000e+00, 1.4141220000e+00, 6.5332900000e-01,
                2.8669200000e-01, 1.1501400000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{3.7900000000e-04, 3.5730000000e-03, 1.9229000000e-02,
                6.9903000000e-02, 1.9938000000e-01, 3.4884800000e-01,
                5.1394000000e-02, -5.5967900000e-01, -3.8458100000e-01,
                4.3892500000e-01, 4.9049800000e-01},
      doubles_t{5.2788600000e+02, 1.5873300000e+02, 6.1555600000e+01,
                2.6883400000e+01, 1.2517500000e+01, 6.0697700000e+00,
                2.9583190000e+00, 1.4141220000e+00, 6.5332900000e-01,
                2.8669200000e-01, 1.1501400000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1501400000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7344000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3909000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5747000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0075000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7760000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3212000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7909000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6768100000e+00}));
    return abs_t(name, 29, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_q_plus_d_cparen_z_29

} // namespace chemcache
