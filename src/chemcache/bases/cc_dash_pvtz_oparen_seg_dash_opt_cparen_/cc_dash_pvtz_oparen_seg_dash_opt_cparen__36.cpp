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

#include "cc_dash_pvtz_oparen_seg_dash_opt_cparen_.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvtz_oparen_seg_dash_opt_cparen__36() {
    // Basis Set name and origin point
    std::string name("cc-pvtz(seg-opt)");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.6000000000e-06, 4.3800000000e-05, 2.3050000000e-04,
                9.7330000000e-04, 3.5335300000e-03, 1.1415690000e-02,
                3.3123220000e-02, 8.5393490000e-02, 1.8662470000e-01,
                3.1388400000e-01, 3.3149949000e-01, 1.5917405000e-01,
                2.7706600000e-02, 1.6032840000e-02, 7.7607300000e-03},
      doubles_t{1.1718113000e+07, 1.7546044000e+06, 3.9928132000e+05,
                1.1308457000e+05, 3.6885925000e+04, 1.3312209000e+04,
                5.1899883000e+03, 2.1516597000e+03, 9.3803251000e+02,
                4.2655732000e+02, 2.0106660000e+02, 9.7097605000e+01,
                4.2998724000e+01, 2.1177075000e+01, 1.0426752000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.0000000000e-08, -2.0000000000e-08, -6.0000000000e-08,
                -5.6400000000e-06, -3.4130000000e-05, -3.1630000000e-04,
                -1.8088900000e-03, -9.8455500000e-03, -3.7389630000e-02,
                -9.7026070000e-02, -5.8426150000e-02, 3.3749687000e-01,
                5.9450348000e-01, 2.4272147000e-01, 1.6695440000e-02},
      doubles_t{1.1718113000e+07, 1.7546044000e+06, 3.9928132000e+05,
                3.6885925000e+04, 1.3312209000e+04, 5.1899883000e+03,
                2.1516597000e+03, 9.3803251000e+02, 4.2655732000e+02,
                2.0106660000e+02, 9.7097605000e+01, 4.2998724000e+01,
                2.1177075000e+01, 1.0426752000e+01, 4.5850080000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.0000000000e-08, -4.0000000000e-08, 8.0000000000e-08,
                -4.5500000000e-06, -3.1800000000e-05, -2.3408000000e-04,
                -1.1358700000e-03, -4.3999200000e-03, -9.3039600000e-03,
                -3.2956100000e-03, 8.0459000000e-04, -1.7106401000e-01,
                -4.6163620000e-02, 6.8590529000e-01, 6.0868959000e-01},
      doubles_t{1.1718113000e+07, 1.7546044000e+06, 3.9928132000e+05,
                1.3312209000e+04, 5.1899883000e+03, 2.1516597000e+03,
                9.3803251000e+02, 4.2655732000e+02, 2.0106660000e+02,
                9.7097605000e+01, 4.2998724000e+01, 2.1177075000e+01,
                1.0426752000e+01, 4.5850080000e+00, 2.1176030000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0705700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{3.0000000000e-08, -6.0000000000e-08, -1.0000000000e-08,
                2.5300000000e-06, 1.5840000000e-05, 9.1740000000e-05,
                3.0496000000e-04, 5.7276000000e-04, -2.7796000000e-04,
                4.5787100000e-03, 4.1596870000e-02, 1.6785240000e-02,
                -1.9636462000e-01, -4.4405111000e-01, 8.7819746000e-01},
      doubles_t{1.1718113000e+07, 1.7546044000e+06, 3.9928132000e+05,
                5.1899883000e+03, 2.1516597000e+03, 9.3803251000e+02,
                4.2655732000e+02, 2.0106660000e+02, 9.7097605000e+01,
                4.2998724000e+01, 2.1177075000e+01, 1.0426752000e+01,
                4.5850080000e+00, 2.1176030000e+00, 3.4922500000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4482100000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.2309000000e-04, 3.6742400000e-03, 1.9931240000e-02,
                7.7426380000e-02, 2.1406906000e-01, 3.8496619000e-01,
                3.6268713000e-01, 1.1606199000e-01, 2.6779500000e-03},
      doubles_t{9.3663090000e+03, 2.2195543000e+03, 7.1945288000e+02,
                2.7346446000e+02, 1.1475225000e+02, 5.1155569000e+01,
                2.3682676000e+01, 1.0875484000e+01, 4.9551310000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.0200000000e-06, -7.7600000000e-05, -1.2105300000e-03,
                -7.5749500000e-03, -2.6225780000e-02, -1.3085090000e-02,
                2.2917328000e-01, 5.5180439000e-01, 3.9629549000e-01},
      doubles_t{9.3663090000e+03, 7.1945288000e+02, 2.7346446000e+02,
                1.1475225000e+02, 5.1155569000e+01, 2.3682676000e+01,
                1.0875484000e+01, 4.9551310000e+00, 2.2172670000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0641000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-6.4000000000e-07, 3.0767000000e-04, 1.6416500000e-03,
                7.0028100000e-03, 1.4845900000e-03, -6.2814980000e-02,
                -1.7880706000e-01, -1.3053072000e-01, 7.4352850000e-01},
      doubles_t{9.3663090000e+03, 2.7346446000e+02, 1.1475225000e+02,
                5.1155569000e+01, 2.3682676000e+01, 1.0875484000e+01,
                4.9551310000e+00, 2.2172670000e+00, 3.2215400000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1761900000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{1.4044000000e-03, 1.2171500000e-02, 5.6391900000e-02,
                           1.6764300000e-01, 3.1773680000e-01, 3.8726470000e-01,
                           2.7280060000e-01},
                 doubles_t{4.4616133000e+02, 1.3396477000e+02, 5.1345907000e+01,
                           2.1916906000e+01, 9.8937250000e+00, 4.4925270000e+00,
                           2.0022930000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0840900000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0060000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6220000000e-01}));
    return abs_t(name, 36, r0, shells.begin(), shells.end());
} // cc_dash_pvtz_oparen_seg_dash_opt_cparen__36

} // namespace chemcache
