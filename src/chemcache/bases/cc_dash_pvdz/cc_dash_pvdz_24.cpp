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

#include "cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvdz_24() {
    // Basis Set name and origin point
    std::string name("cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.1286670000e-06, 3.2107670000e-05,  1.6884160000e-04,
                7.1285200000e-04, 2.5893250000e-03,  8.3773500000e-03,
                2.4417250000e-02, 6.3651350000e-02,  1.4276180000e-01,
                2.5412750000e-01, 3.0095120000e-01,  1.7665130000e-01,
                6.9367090000e-02, 1.1795790000e-01,  8.9161870000e-02,
                1.1036300000e-02, -3.5460480000e-04, 1.0573110000e-04,
                1.1146400000e-05, 2.6613870000e-05},
      doubles_t{6.1771940000e+06, 9.2492950000e+05, 2.1048650000e+05,
                5.9620050000e+04, 1.9450760000e+04, 7.0220560000e+03,
                2.7387630000e+03, 1.1358140000e+03, 4.9509230000e+02,
                2.2474870000e+02, 1.0538360000e+02, 5.0193590000e+01,
                2.2249570000e+01, 1.0982650000e+01, 5.3836650000e+00,
                2.3436850000e+00, 1.1052020000e+00, 4.8784800000e-01,
                8.9599000000e-02, 3.3423000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.3017720000e-06, -1.7895360000e-05, -9.4161740000e-05,
                -3.9750740000e-04, -1.4470250000e-03, -4.6946220000e-03,
                -1.3823870000e-02, -3.6746430000e-02, -8.6471850000e-02,
                -1.6967350000e-01, -2.5070890000e-01, -1.9611560000e-01,
                1.4572440000e-01,  5.4667060000e-01,  3.9794340000e-01,
                5.2770070000e-02,  -4.3745370000e-03, 3.2040350000e-04,
                -5.1420770000e-05, 1.5841340000e-04},
      doubles_t{6.1771940000e+06, 9.2492950000e+05, 2.1048650000e+05,
                5.9620050000e+04, 1.9450760000e+04, 7.0220560000e+03,
                2.7387630000e+03, 1.1358140000e+03, 4.9509230000e+02,
                2.2474870000e+02, 1.0538360000e+02, 5.0193590000e+01,
                2.2249570000e+01, 1.0982650000e+01, 5.3836650000e+00,
                2.3436850000e+00, 1.1052020000e+00, 4.8784800000e-01,
                8.9599000000e-02, 3.3423000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.8629570000e-07,  3.7766450000e-06,  1.9906640000e-05,
                8.3891640000e-05,  3.0657060000e-04,  9.9441070000e-04,
                2.9619590000e-03,  7.9694730000e-03,  1.9550170000e-02,
                4.0850350000e-02,  6.9290030000e-02,  6.1469840000e-02,
                -6.9813020000e-02, -3.5175970000e-01, -3.8286290000e-01,
                2.6764010000e-01,  7.1759500000e-01,  3.0208140000e-01,
                7.7495140000e-03,  2.6960960000e-04},
      doubles_t{6.1771940000e+06, 9.2492950000e+05, 2.1048650000e+05,
                5.9620050000e+04, 1.9450760000e+04, 7.0220560000e+03,
                2.7387630000e+03, 1.1358140000e+03, 4.9509230000e+02,
                2.2474870000e+02, 1.0538360000e+02, 5.0193590000e+01,
                2.2249570000e+01, 1.0982650000e+01, 5.3836650000e+00,
                2.3436850000e+00, 1.1052020000e+00, 4.8784800000e-01,
                8.9599000000e-02, 3.3423000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1024510000e-07, -8.5302330000e-07, -4.5203580000e-06,
                -1.8916120000e-05, -6.9743440000e-05, -2.2378670000e-04,
                -6.7545030000e-04, -1.7893460000e-03, -4.4778580000e-03,
                -9.1401440000e-03, -1.6105620000e-02, -1.3348700000e-02,
                1.4260270000e-02,  8.9316900000e-02,  8.8852790000e-02,
                -6.3687760000e-02, -2.7832620000e-01, -1.8300710000e-01,
                6.7909370000e-01,  4.6729530000e-01},
      doubles_t{6.1771940000e+06, 9.2492950000e+05, 2.1048650000e+05,
                5.9620050000e+04, 1.9450760000e+04, 7.0220560000e+03,
                2.7387630000e+03, 1.1358140000e+03, 4.9509230000e+02,
                2.2474870000e+02, 1.0538360000e+02, 5.0193590000e+01,
                2.2249570000e+01, 1.0982650000e+01, 5.3836650000e+00,
                2.3436850000e+00, 1.1052020000e+00, 4.8784800000e-01,
                8.9599000000e-02, 3.3423000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3423000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.1798930000e-07,  1.6129400000e-06,  9.1118420000e-06,
                3.5006450000e-05,  1.4353150000e-04,  4.0358960000e-04,
                1.4251770000e-03,  3.1140090000e-03,  9.8144490000e-03,
                1.4746980000e-02,  3.9115120000e-02,  9.1708880000e-03,
                1.5598780000e-02,  -2.8168440000e-01, -6.8952610000e-03,
                -1.7697810000e-01, 1.4430610000e+00,  -1.0293180000e+00,
                -1.3076670000e+00, 1.5038420000e+00},
      doubles_t{6.1771940000e+06, 9.2492950000e+05, 2.1048650000e+05,
                5.9620050000e+04, 1.9450760000e+04, 7.0220560000e+03,
                2.7387630000e+03, 1.1358140000e+03, 4.9509230000e+02,
                2.2474870000e+02, 1.0538360000e+02, 5.0193590000e+01,
                2.2249570000e+01, 1.0982650000e+01, 5.3836650000e+00,
                2.3436850000e+00, 1.1052020000e+00, 4.8784800000e-01,
                8.9599000000e-02, 3.3423000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.4000000000e-05, 3.8900000000e-04, 2.2410000000e-03,
                9.8210000000e-03, 3.4471000000e-02, 9.7460000000e-02,
                2.1198500000e-01, 3.3399000000e-01, 3.3013700000e-01,
                1.5222700000e-01, 2.0425000000e-02, -1.3600000000e-03,
                -1.1950000000e-03, -1.9700000000e-04, 2.3000000000e-05,
                -9.0000000000e-06},
      doubles_t{1.4454200000e+04, 3.4216760000e+03, 1.1113870000e+03,
                4.2519180000e+02, 1.8026230000e+02, 8.2061170000e+01,
                3.9297260000e+01, 1.9419590000e+01, 9.8288990000e+00,
                5.0168100000e+00, 2.4870910000e+00, 1.1987800000e+00,
                5.5869500000e-01, 2.0892400000e-01, 8.4608000000e-02,
                3.3258000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5000000000e-05, -1.3500000000e-04, -7.7700000000e-04,
                -3.4270000000e-03, -1.2189000000e-02, -3.5388000000e-02,
                -7.9915000000e-02, -1.3233500000e-01, -1.3540100000e-01,
                3.2008000000e-02, 3.3384900000e-01, 4.6177300000e-01,
                2.8129000000e-01, 4.1843000000e-02, -4.0020000000e-03,
                1.5210000000e-03},
      doubles_t{1.4454200000e+04, 3.4216760000e+03, 1.1113870000e+03,
                4.2519180000e+02, 1.8026230000e+02, 8.2061170000e+01,
                3.9297260000e+01, 1.9419590000e+01, 9.8288990000e+00,
                5.0168100000e+00, 2.4870910000e+00, 1.1987800000e+00,
                5.5869500000e-01, 2.0892400000e-01, 8.4608000000e-02,
                3.3258000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0000000000e-06, 4.0000000000e-05, 2.2900000000e-04,
                1.0190000000e-03, 3.6020000000e-03, 1.0550000000e-02,
                2.3702000000e-02, 3.9988000000e-02, 4.0437000000e-02,
                -1.2074000000e-02, -1.1893900000e-01, -1.7810000000e-01,
                -1.2386500000e-01, 4.2972200000e-01, 6.5078600000e-01,
                6.4171000000e-02},
      doubles_t{1.4454200000e+04, 3.4216760000e+03, 1.1113870000e+03,
                4.2519180000e+02, 1.8026230000e+02, 8.2061170000e+01,
                3.9297260000e+01, 1.9419590000e+01, 9.8288990000e+00,
                5.0168100000e+00, 2.4870910000e+00, 1.1987800000e+00,
                5.5869500000e-01, 2.0892400000e-01, 8.4608000000e-02,
                3.3258000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0000000000e-06, 3.2000000000e-05, 1.8500000000e-04,
                8.1000000000e-04, 2.9060000000e-03, 8.3910000000e-03,
                1.9193000000e-02, 3.1564000000e-02, 3.3417000000e-02,
                -1.2907000000e-02, -9.3659000000e-02, -1.4997700000e-01,
                -6.7234000000e-02, 2.7075900000e-01, 5.7580700000e-01,
                3.0112100000e-01},
      doubles_t{1.4454200000e+04, 3.4216760000e+03, 1.1113870000e+03,
                4.2519180000e+02, 1.8026230000e+02, 8.2061170000e+01,
                3.9297260000e+01, 1.9419590000e+01, 9.8288990000e+00,
                5.0168100000e+00, 2.4870910000e+00, 1.1987800000e+00,
                5.5869500000e-01, 2.0892400000e-01, 8.4608000000e-02,
                3.3258000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3258000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.6210000000e-03, 2.5766000000e-02, 9.7556000000e-02,
                           2.3631200000e-01, 3.5828600000e-01, 3.6854300000e-01,
                           2.3549400000e-01, 5.3156000000e-02},
                 doubles_t{8.8576800000e+01, 2.6204500000e+01, 9.5174700000e+00,
                           3.8224800000e+00, 1.5751200000e+00, 6.2892800000e-01,
                           2.3442400000e-01, 7.6815000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-4.1220000000e-03, -2.9307000000e-02, -1.1506200000e-01,
                -2.7306800000e-01, -3.1442300000e-01, 4.2097000000e-02,
                5.9140300000e-01, 3.5821500000e-01},
      doubles_t{8.8576800000e+01, 2.6204500000e+01, 9.5174700000e+00,
                3.8224800000e+00, 1.5751200000e+00, 6.2892800000e-01,
                2.3442400000e-01, 7.6815000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6815000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3, doubles_t{4.2354500000e-01, 7.7411400000e-01},
                 doubles_t{2.2211000000e+00, 5.2310000000e-01}));
    return abs_t(name, 24, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_24

} // namespace chemcache
