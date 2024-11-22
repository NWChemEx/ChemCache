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

#include "aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_34() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(d+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0040000000e-04, 1.5554000000e-03, 8.0872000000e-03,
                3.3034400000e-02, 1.0829240000e-01, 2.7033610000e-01,
                4.2062360000e-01, 2.8159220000e-01, 3.0811000000e-02,
                -7.1617000000e-03, 2.6022000000e-03, -1.2583000000e-03,
                3.4650000000e-04, -1.5030000000e-04},
      doubles_t{5.9899000000e+05, 8.9783000000e+04, 2.0435000000e+04,
                5.7869000000e+03, 1.8873000000e+03, 6.8097000000e+02,
                2.6539000000e+02, 1.0863000000e+02, 3.3760000000e+01,
                1.4465000000e+01, 4.3890000000e+00, 1.8783000000e+00,
                3.5859000000e-01, 1.3649000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.2900000000e-05, -4.8500000000e-04, -2.5644000000e-03,
                -1.0476100000e-02, -3.6722300000e-02, -9.9922500000e-02,
                -2.1419730000e-01, -1.8365930000e-01, 4.6754540000e-01,
                6.4147400000e-01, 7.9256900000e-02, -1.4269700000e-02,
                3.3792000000e-03, -1.4537000000e-03},
      doubles_t{5.9899000000e+05, 8.9783000000e+04, 2.0435000000e+04,
                5.7869000000e+03, 1.8873000000e+03, 6.8097000000e+02,
                2.6539000000e+02, 1.0863000000e+02, 3.3760000000e+01,
                1.4465000000e+01, 4.3890000000e+00, 1.8783000000e+00,
                3.5859000000e-01, 1.3649000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4700000000e-05, 1.9130000000e-04, 1.0068000000e-03,
                4.1514000000e-03, 1.4499100000e-02, 4.0565800000e-02,
                8.8536400000e-02, 8.5421200000e-02, -2.9325810000e-01,
                -5.5707270000e-01, 5.2614360000e-01, 7.3203710000e-01,
                3.8824600000e-02, -1.0503600000e-02},
      doubles_t{5.9899000000e+05, 8.9783000000e+04, 2.0435000000e+04,
                5.7869000000e+03, 1.8873000000e+03, 6.8097000000e+02,
                2.6539000000e+02, 1.0863000000e+02, 3.3760000000e+01,
                1.4465000000e+01, 4.3890000000e+00, 1.8783000000e+00,
                3.5859000000e-01, 1.3649000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.2000000000e-06, -5.5900000000e-05, -2.9380000000e-04,
                -1.2136000000e-03, -4.2340000000e-03, -1.1903500000e-02,
                -2.6020600000e-02, -2.5614800000e-02, 9.1942700000e-02,
                1.8387000000e-01, -2.1884610000e-01, -4.8965240000e-01,
                6.7758180000e-01, 5.2967210000e-01},
      doubles_t{5.9899000000e+05, 8.9783000000e+04, 2.0435000000e+04,
                5.7869000000e+03, 1.8873000000e+03, 6.8097000000e+02,
                2.6539000000e+02, 1.0863000000e+02, 3.3760000000e+01,
                1.4465000000e+01, 4.3890000000e+00, 1.8783000000e+00,
                3.5859000000e-01, 1.3649000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3649000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8747000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.4127000000e-03, 1.1858800000e-02, 5.9515300000e-02,
                           1.9722010000e-01, 4.0074390000e-01, 3.9947400000e-01,
                           1.1533640000e-01, 2.2190000000e-04, 2.2838000000e-03,
                           -4.7560000000e-04, 1.5160000000e-04},
                 doubles_t{4.1356000000e+03, 9.8034000000e+02, 3.1635000000e+02,
                           1.1925000000e+02, 4.9068000000e+01, 2.1212000000e+01,
                           8.9462000000e+00, 3.8236000000e+00, 1.5883000000e+00,
                           4.0969000000e-01, 1.2459000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-5.6100000000e-04, -4.7340000000e-03, -2.4350400000e-02,
                -8.4107100000e-02, -1.8413840000e-01, -1.7350040000e-01,
                2.1672630000e-01, 5.8500990000e-01, 3.4168160000e-01,
                1.9912500000e-02, -2.6131000000e-03},
      doubles_t{4.1356000000e+03, 9.8034000000e+02, 3.1635000000e+02,
                1.1925000000e+02, 4.9068000000e+01, 2.1212000000e+01,
                8.9462000000e+00, 3.8236000000e+00, 1.5883000000e+00,
                4.0969000000e-01, 1.2459000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.3660000000e-04, 1.1308000000e-03, 5.9581000000e-03,
                2.0186600000e-02, 4.6193900000e-02, 3.9405000000e-02,
                -5.9284600000e-02, -2.0146630000e-01, -6.8782100000e-02,
                5.5959440000e-01, 5.7097840000e-01},
      doubles_t{4.1356000000e+03, 9.8034000000e+02, 3.1635000000e+02,
                1.1925000000e+02, 4.9068000000e+01, 2.1212000000e+01,
                8.9462000000e+00, 3.8236000000e+00, 1.5883000000e+00,
                4.0969000000e-01, 1.2459000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2459000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5492000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.3498200000e-02, 1.3751830000e-01, 3.6648240000e-01,
                           4.8747170000e-01, 2.7657690000e-01},
                 doubles_t{9.4472000000e+01, 2.7180000000e+01, 9.5068000000e+00,
                           3.4168000000e+00, 1.1479000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6820000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2830000000e-01}));
    return abs_t(name, 34, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_d_plus_d_cparen_z_34

} // namespace chemcache
