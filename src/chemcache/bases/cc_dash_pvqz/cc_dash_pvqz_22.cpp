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

abs_t cc_dash_pvqz_22() {
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
      doubles_t{1.9678840000e-06, 1.5301980000e-05,  8.0504030000e-05,
                3.4003020000e-04, 1.2383490000e-03,  4.0284750000e-03,
                1.1925080000e-02, 3.2156010000e-02,  7.7786200000e-02,
                1.6120300000e-01, 2.6188260000e-01,  2.7665410000e-01,
                1.4429680000e-01, 6.9931310000e-02,  1.2594490000e-01,
                9.0258300000e-02, 1.2004010000e-02,  -3.9761680000e-04,
                1.4254440000e-04, -4.1966640000e-05, 4.8242130000e-05,
                -4.6768980000e-06},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1300070000e-06, -8.7894160000e-06, -4.6226230000e-05,
                -1.9540470000e-04, -7.1170000000e-04, -2.3217550000e-03,
                -6.8995110000e-03, -1.8837650000e-02, -4.6639760000e-02,
                -1.0207640000e-01, -1.8486930000e-01, -2.4796010000e-01,
                -1.6348630000e-01, 1.7791470000e-01,  5.3517380000e-01,
                3.7528220000e-01,  5.2170380000e-02,  -3.8631920000e-03,
                2.8453740000e-04,  -3.5004140000e-04, 3.1929060000e-04,
                -5.6139700000e-05},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.2527350000e-07,  1.7536310000e-06,  9.2130520000e-06,
                3.9012420000e-05,  1.4191840000e-04,  4.6486350000e-04,
                1.3835440000e-03,  3.8273540000e-03,  9.6391910000e-03,
                2.2117190000e-02,  4.3149540000e-02,  6.7034810000e-02,
                4.9464240000e-02,  -8.3334560000e-02, -3.3844350000e-01,
                -3.6740690000e-01, 2.2470560000e-01,  6.8788150000e-01,
                3.5944350000e-01,  1.9336950000e-02,  -5.9978190000e-03,
                2.8280640000e-03},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.3184330000e-08, -4.1398550000e-07, -2.1751530000e-06,
                -9.2095670000e-06, -3.3507610000e-05, -1.0973930000e-04,
                -3.2670060000e-04, -9.0365670000e-04, -2.2772570000e-03,
                -5.2267560000e-03, -1.0217070000e-02, -1.5914680000e-02,
                -1.1828140000e-02, 2.0119990000e-02,  8.5265840000e-02,
                9.6934890000e-02,  -7.0826720000e-02, -2.3597400000e-01,
                -2.6296130000e-01, 3.0112140000e-01,  6.1250740000e-01,
                2.5713680000e-01},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.7118720000e-08, -7.7339700000e-07, -3.9310520000e-06,
                -1.7383230000e-05, -5.9865360000e-05, -2.0958530000e-04,
                -5.7564820000e-04, -1.7515950000e-03, -3.9351200000e-03,
                -1.0384900000e-02, -1.6911000000e-02, -3.4109890000e-02,
                -1.2204370000e-02, 1.3975350000e-02,  2.1762120000e-01,
                9.5327080000e-02,  1.7825240000e-02,  -1.0163420000e+00,
                3.1070260000e-01,  1.9917860000e+00,  -1.1429550000e+00,
                -6.2069230000e-01},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.5336440000e-07, -1.2342950000e-06, -6.1771880000e-06,
                -2.7872530000e-05, -9.3551990000e-05, -3.3783520000e-04,
                -8.9354750000e-04, -2.8423530000e-03, -6.0523320000e-03,
                -1.7046900000e-02, -2.5498060000e-02, -5.7911490000e-02,
                -1.2717910000e-02, 6.2607620000e-03,  4.0803060000e-01,
                1.2842580000e-01,  -1.2312860000e-01, -2.2673310000e+00,
                2.8075070000e+00,  -3.0280090000e-01, -2.0712270000e+00,
                1.8036250000e+00},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7645260000e-07, -1.3144010000e-06, -7.3554640000e-06,
                -2.8637620000e-05, -1.1565510000e-04, -3.3307320000e-04,
                -1.1551560000e-03, -2.6598320000e-03, -8.3320360000e-03,
                -1.4626110000e-02, -4.0308900000e-02, -3.7458420000e-02,
                -7.3759590000e-02, 1.4432680000e-01,  2.1190310000e-01,
                7.4741870000e-01,  -2.0310800000e+00, 1.9226930000e-01,
                2.2877680000e+00,  -5.1681350000e+00, 6.4603230000e+00,
                -2.9155510000e+00},
      doubles_t{
        9.3170340000e+06, 1.3949710000e+06, 3.1745310000e+05, 8.9921000000e+04,
        2.9338160000e+04, 1.0592660000e+04, 4.1320560000e+03, 1.7142350000e+03,
        7.4791510000e+02, 3.4019910000e+02, 1.6020600000e+02, 7.7588510000e+01,
        3.8044890000e+01, 1.7353620000e+01, 8.6420620000e+00, 4.3044040000e+00,
        1.9377180000e+00, 9.3659400000e-01, 4.3262400000e-01, 1.0735200000e-01,
        5.2034000000e-02, 2.4038000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.2000000000e-05, 1.0400000000e-04, 6.0700000000e-04,
                2.7510000000e-03, 1.0266000000e-02, 3.2180000000e-02,
                8.4486000000e-02, 1.7833200000e-01, 2.8857400000e-01,
                3.2875800000e-01, 2.1154100000e-01, 5.5519000000e-02,
                2.0200000000e-03, -1.4290000000e-03, -9.4100000000e-04,
                -5.3000000000e-05, -1.7000000000e-05, 3.0000000000e-06},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-4.0000000000e-06, -3.5000000000e-05, -2.0300000000e-04,
                -9.2500000000e-04, -3.4700000000e-03, -1.1020000000e-02,
                -2.9577000000e-02, -6.4574000000e-02, -1.0845800000e-01,
                -1.3191000000e-01, -5.3089000000e-02, 1.7812500000e-01,
                3.9570100000e-01, 3.9272800000e-01, 1.7447700000e-01,
                1.8946000000e-02, -3.2800000000e-04, 3.7600000000e-04},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 1.0000000000e-05, 6.1000000000e-05,
                2.7500000000e-04, 1.0340000000e-03, 3.2830000000e-03,
                8.8290000000e-03, 1.9308000000e-02, 3.2612000000e-02,
                3.9977000000e-02, 1.5299000000e-02, -6.2620000000e-02,
                -1.3828800000e-01, -1.7634800000e-01, -5.5404000000e-02,
                5.5929100000e-01, 5.3311500000e-01, 3.8172000000e-02},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0000000000e-06, 2.0000000000e-05, 1.2600000000e-04,
                5.2600000000e-04, 2.1510000000e-03, 6.2570000000e-03,
                1.8455000000e-02, 3.6247000000e-02, 7.0920000000e-02,
                7.1827000000e-02, 5.9201000000e-02, -1.8427200000e-01,
                -3.2965300000e-01, -5.4576400000e-01, 1.1454310000e+00,
                5.5262500000e-01, -9.6426000000e-01, -8.8240000000e-02},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.0000000000e-06, -2.8000000000e-05, -1.3400000000e-04,
                -7.3600000000e-04, -2.2930000000e-03, -8.7520000000e-03,
                -1.9226000000e-02, -5.2589000000e-02, -6.5894000000e-02,
                -1.4811200000e-01, -2.4124000000e-02, 1.7697900000e-01,
                1.1584280000e+00, -9.6345800000e-01, -1.3245080000e+00,
                2.2913730000e+00, -1.1129320000e+00, -2.1563100000e-01},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0000000000e-06, 8.0000000000e-06, 4.9000000000e-05,
                2.2100000000e-04, 8.3600000000e-04, 2.6410000000e-03,
                7.1510000000e-03, 1.5530000000e-02, 2.6479000000e-02,
                3.1980000000e-02, 1.2561000000e-02, -5.2994000000e-02,
                -1.1303900000e-01, -1.4369800000e-01, -1.0336000000e-02,
                3.5473600000e-01, 5.4810100000e-01, 2.2529600000e-01},
      doubles_t{2.5537330000e+04, 6.0422430000e+03, 1.9626530000e+03,
                7.5173380000e+02, 3.1988330000e+02, 1.4643000000e+02,
                7.0748610000e+01, 3.5590880000e+01, 1.8368930000e+01,
                9.6586550000e+00, 5.1367300000e+00, 2.6988720000e+00,
                1.3897030000e+00, 7.0316300000e-01, 3.4412300000e-01,
                1.4136100000e-01, 6.0863000000e-02, 2.5907000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5907000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.8200000000e-04, 3.3680000000e-03, 1.6242000000e-02,
                           5.2377000000e-02, 1.2992800000e-01, 2.3123300000e-01,
                           3.0356100000e-01, 3.0926800000e-01, 2.3164400000e-01,
                           1.0154000000e-01, 1.1684000000e-02},
                 doubles_t{1.9694900000e+02, 5.8837200000e+01, 2.2565400000e+01,
                           9.5756600000e+00, 4.3185900000e+00, 2.0232700000e+00,
                           9.4677800000e-01, 4.3387600000e-01, 1.9227400000e-01,
                           8.1502000000e-02, 3.2285000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.8800000000e-04, -3.4280000000e-03, -1.6598000000e-02,
                -5.4172000000e-02, -1.3666000000e-01, -2.2944900000e-01,
                -2.5499500000e-01, -4.7661000000e-02, 3.9488000000e-01,
                4.9645300000e-01, 1.2575600000e-01},
      doubles_t{1.9694900000e+02, 5.8837200000e+01, 2.2565400000e+01,
                9.5756600000e+00, 4.3185900000e+00, 2.0232700000e+00,
                9.4677800000e-01, 4.3387600000e-01, 1.9227400000e-01,
                8.1502000000e-02, 3.2285000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-6.2200000000e-04, -5.4930000000e-03, -2.6968000000e-02,
                -8.9533000000e-02, -2.2837400000e-01, -3.3125800000e-01,
                -1.0634900000e-01, 6.3103600000e-01, 3.9662100000e-01,
                -6.1800900000e-01, -2.6959200000e-01},
      doubles_t{1.9694900000e+02, 5.8837200000e+01, 2.2565400000e+01,
                9.5756600000e+00, 4.3185900000e+00, 2.0232700000e+00,
                9.4677800000e-01, 4.3387600000e-01, 1.9227400000e-01,
                8.1502000000e-02, 3.2285000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-6.9600000000e-04, -7.0890000000e-03, -3.1013000000e-02,
                -1.2005500000e-01, -2.7891000000e-01, -3.7896100000e-01,
                5.7374300000e-01, 7.2795000000e-01, -1.2976560000e+00,
                3.2037600000e-01, 5.2301000000e-01},
      doubles_t{1.9694900000e+02, 5.8837200000e+01, 2.2565400000e+01,
                9.5756600000e+00, 4.3185900000e+00, 2.0232700000e+00,
                9.4677800000e-01, 4.3387600000e-01, 1.9227400000e-01,
                8.1502000000e-02, 3.2285000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2285000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2626000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4180000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8210000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5732000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0020000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7210000000e-01}));
    return abs_t(name, 22, r0, shells.begin(), shells.end());
} // cc_dash_pvqz_22

} // namespace chemcache
