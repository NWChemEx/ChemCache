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

#include "cc_dash_pvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvtz_30() {
    // Basis Set name and origin point
    std::string name("cc-pvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{8.5492410000e-06,  6.6474100000e-05,  3.4949620000e-04,
                1.4738320000e-03,  5.3383300000e-03,  1.7127080000e-02,
                4.8940850000e-02,  1.2179340000e-01,  2.4765890000e-01,
                3.5824310000e-01,  2.7981740000e-01,  6.8574910000e-02,
                -1.3110920000e-03, 1.9140010000e-03,  -8.7592200000e-04,
                3.7400960000e-04,  -1.4013990000e-04, 4.7571320000e-05,
                -3.6427110000e-05, 1.1532480000e-05},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.6400690000e-06, -2.0527200000e-05, -1.0798590000e-04,
                -4.5585770000e-04, -1.6577580000e-03, -5.3684920000e-03,
                -1.5712490000e-02, -4.1225580000e-02, -9.4064590000e-02,
                -1.7199540000e-01, -1.9585230000e-01, 4.5329070000e-02,
                5.2444420000e-01,  5.0061420000e-01,  8.9455270000e-02,
                -2.1462620000e-03, 2.1121130000e-03,  -4.1339800000e-04,
                3.2097520000e-04,  -1.0161400000e-04},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.9671030000e-07,  7.7541630000e-06,  4.0760190000e-05,
                1.7228110000e-04,  6.2593700000e-04,  2.0328550000e-03,
                5.9546460000e-03,  1.5766400000e-02,  3.6376380000e-02,
                6.8923430000e-02,  8.2380930000e-02,  -2.0113600000e-02,
                -3.2525260000e-01, -4.6028990000e-01, 1.6355460000e-01,
                7.2971180000e-01,  3.7697510000e-01,  1.4332240000e-02,
                -6.6712100000e-03, 1.7662140000e-03},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9958180000e-07,  1.5529730000e-06,  8.1612590000e-06,
                3.4507470000e-05,  1.2532750000e-04,  4.0729900000e-04,
                1.1927340000e-03,  3.1631400000e-03,  7.3039420000e-03,
                1.3912790000e-02,  1.6706200000e-02,  -4.0355860000e-03,
                -6.9688610000e-02, -1.0301050000e-01, 4.4714420000e-02,
                2.1500270000e-01,  2.2201630000e-01,  -3.1147760000e-01,
                -5.6934290000e-01, -2.6784400000e-01},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{8.0407910000e-07,  6.2829360000e-06,  3.2818680000e-05,
                1.3988580000e-04,  5.0298360000e-04,  1.6554190000e-03,
                4.7786770000e-03,  1.2924790000e-02,  2.9258310000e-02,
                5.7918160000e-02,  6.6406810000e-02,  -7.3889660000e-03,
                -3.3299890000e-01, -5.9178650000e-01, 9.0114060000e-01,
                1.5859510000e+00,  -2.7880080000e+00, 2.0718840000e+00,
                -6.0120250000e-01, -5.9371840000e-01},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.4359100000e-07, -4.3368940000e-06, -2.1975720000e-05,
                -9.7473920000e-05, -3.3316150000e-04, -1.1661920000e-03,
                -3.1193080000e-03, -9.2395040000e-03, -1.8554710000e-02,
                -4.2811890000e-02, -3.5710950000e-02, -1.6383500000e-02,
                2.6446640000e-01,  2.0865880000e-01,  -1.7743820000e-02,
                -1.3538730000e+00, 8.1829260000e-01,  1.6950360000e+00,
                -1.3886560000e+00, -2.1889000000e-01},
      doubles_t{5.8200210000e+06, 8.7152340000e+05, 1.9833500000e+05,
                5.6176310000e+04, 1.8325820000e+04, 6.6149550000e+03,
                2.5791990000e+03, 1.0688490000e+03, 4.6510450000e+02,
                2.1041300000e+02, 9.7616290000e+01, 4.4380200000e+01,
                2.1423080000e+01, 1.0308910000e+01, 4.5536450000e+00,
                2.1328210000e+00, 9.2969700000e-01, 1.9214700000e-01,
                8.7595000000e-02, 3.7702000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7702000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.1000000000e-05, 3.6100000000e-04, 2.0880000000e-03,
                9.2210000000e-03, 3.2773000000e-02, 9.4179000000e-02,
                2.0913200000e-01, 3.3456900000e-01, 3.3035900000e-01,
                1.5234700000e-01, 2.2984000000e-02, 1.6070000000e-03,
                4.6800000000e-04, 6.6000000000e-05, -2.0000000000e-06},
      doubles_t{2.4411980000e+04, 5.7785180000e+03, 1.8768620000e+03,
                7.1823610000e+02, 3.0483270000e+02, 1.3904530000e+02,
                6.6804170000e+01, 3.3206990000e+01, 1.6928160000e+01,
                8.6962290000e+00, 4.3505100000e+00, 2.1165230000e+00,
                9.9538700000e-01, 3.7811200000e-01, 1.3457900000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5000000000e-05, -1.3500000000e-04, -7.8200000000e-04,
                -3.4780000000e-03, -1.2520000000e-02, -3.7016000000e-02,
                -8.5559000000e-02, -1.4471800000e-01, -1.4344200000e-01,
                4.3595000000e-02, 3.4888800000e-01, 4.5386500000e-01,
                2.6859400000e-01, 3.8868000000e-02, -2.4920000000e-03,
                1.0140000000e-03},
      doubles_t{2.4411980000e+04, 5.7785180000e+03, 1.8768620000e+03,
                7.1823610000e+02, 3.0483270000e+02, 1.3904530000e+02,
                6.6804170000e+01, 3.3206990000e+01, 1.6928160000e+01,
                8.6962290000e+00, 4.3505100000e+00, 2.1165230000e+00,
                9.9538700000e-01, 3.7811200000e-01, 1.3457900000e-01,
                4.6282000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.2000000000e-05, 9.6000000000e-05, 5.9400000000e-04,
                2.4840000000e-03, 9.5370000000e-03, 2.6479000000e-02,
                6.6366000000e-02, 1.0245800000e-01, 1.3868300000e-01,
                -8.0140000000e-02, -4.9606900000e-01, -4.6351000000e-01,
                8.7453100000e-01, 6.2979000000e-01, -8.1168600000e-01,
                -1.0894800000e-01},
      doubles_t{2.4411980000e+04, 5.7785180000e+03, 1.8768620000e+03,
                7.1823610000e+02, 3.0483270000e+02, 1.3904530000e+02,
                6.6804170000e+01, 3.3206990000e+01, 1.6928160000e+01,
                8.6962290000e+00, 4.3505100000e+00, 2.1165230000e+00,
                9.9538700000e-01, 3.7811200000e-01, 1.3457900000e-01,
                4.6282000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.0000000000e-06, 4.2000000000e-05, 2.3800000000e-04,
                1.0880000000e-03, 3.8210000000e-03, 1.1644000000e-02,
                2.6167000000e-02, 4.6750000000e-02, 4.3309000000e-02,
                -1.3429000000e-02, -1.5389700000e-01, -1.6741300000e-01,
                -8.4995000000e-02, 4.5081300000e-01, 6.4086900000e-01,
                5.4172000000e-02},
      doubles_t{2.4411980000e+04, 5.7785180000e+03, 1.8768620000e+03,
                7.1823610000e+02, 3.0483270000e+02, 1.3904530000e+02,
                6.6804170000e+01, 3.3206990000e+01, 1.6928160000e+01,
                8.6962290000e+00, 4.3505100000e+00, 2.1165230000e+00,
                9.9538700000e-01, 3.7811200000e-01, 1.3457900000e-01,
                4.6282000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06, 2.5000000000e-05, 1.4400000000e-04,
                6.4500000000e-04, 2.3110000000e-03, 6.8980000000e-03,
                1.5882000000e-02, 2.7350000000e-02, 2.6621000000e-02,
                -1.0858000000e-02, -7.9853000000e-02, -1.0612700000e-01,
                -6.8883000000e-02, 1.8438500000e-01, 5.6178800000e-01,
                4.1441600000e-01},
      doubles_t{2.4411980000e+04, 5.7785180000e+03, 1.8768620000e+03,
                7.1823610000e+02, 3.0483270000e+02, 1.3904530000e+02,
                6.6804170000e+01, 3.3206990000e+01, 1.6928160000e+01,
                8.6962290000e+00, 4.3505100000e+00, 2.1165230000e+00,
                9.9538700000e-01, 3.7811200000e-01, 1.3457900000e-01,
                4.6282000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6282000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.3420000000e-03, 1.8606000000e-02, 7.7102000000e-02,
                           2.0202600000e-01, 3.2945400000e-01, 3.6097600000e-01,
                           2.7165700000e-01, 1.0498100000e-01},
                 doubles_t{2.0561770000e+02, 6.1449810000e+01, 2.3056890000e+01,
                           9.5777390000e+00, 4.1337340000e+00, 1.7475180000e+00,
                           6.9956000000e-01, 2.5160800000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{3.2790000000e-03, 2.6176000000e-02, 1.1136700000e-01,
                3.0458100000e-01, 3.8629900000e-01, -5.8375000000e-02,
                -5.3887600000e-01, -3.4547300000e-01},
      doubles_t{2.0561770000e+02, 6.1449810000e+01, 2.3056890000e+01,
                9.5777390000e+00, 4.1337340000e+00, 1.7475180000e+00,
                6.9956000000e-01, 2.5160800000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{3.7400000000e-03, 3.1825000000e-02, 1.3222900000e-01,
                4.2455000000e-01, 1.2037000000e-01, -7.6266100000e-01,
                -1.1282300000e-01, 8.0962300000e-01},
      doubles_t{2.0561770000e+02, 6.1449810000e+01, 2.3056890000e+01,
                9.5777390000e+00, 4.1337340000e+00, 1.7475180000e+00,
                6.9956000000e-01, 2.5160800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5160800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7922000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4851000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1144000000e+00}));
    return abs_t(name, 30, r0, shells.begin(), shells.end());
} // cc_dash_pvtz_30

} // namespace chemcache
