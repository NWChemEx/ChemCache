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

#include "aug_dash_cc_dash_pwcvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pwcvdz_16() {
    // Basis Set name and origin point
    std::string name("aug-cc-pwcvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4763500000e-04, 1.9202600000e-03, 9.9619200000e-03,
                4.0297500000e-02, 1.2860400000e-01, 3.0348000000e-01,
                4.2143200000e-01, 2.3078100000e-01, 1.7897100000e-02,
                -2.9751600000e-03, 8.4952200000e-04, -3.6793600000e-04},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                1.0710000000e+03, 3.4980000000e+02, 1.2630000000e+02,
                4.9260000000e+01, 2.0160000000e+01, 5.7200000000e+00,
                2.1820000000e+00, 4.3270000000e-01, 1.5700000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.8703900000e-05, -5.2768100000e-04, -2.7967100000e-03,
                -1.1265100000e-02, -3.8883400000e-02, -9.9502500000e-02,
                -1.9974000000e-01, -1.2336000000e-01, 5.1319400000e-01,
                6.0712000000e-01, 3.9675300000e-02, -9.4686400000e-03},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                1.0710000000e+03, 3.4980000000e+02, 1.2630000000e+02,
                4.9260000000e+01, 2.0160000000e+01, 5.7200000000e+00,
                2.1820000000e+00, 4.3270000000e-01, 1.5700000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9907700000e-05, 1.5348300000e-04, 8.0950300000e-04,
                3.2897400000e-03, 1.1296700000e-02, 2.9638500000e-02,
                5.9985100000e-02, 4.1324800000e-02, -2.0747400000e-01,
                -3.9288900000e-01, 6.3284000000e-01, 5.5692400000e-01},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                1.0710000000e+03, 3.4980000000e+02, 1.2630000000e+02,
                4.9260000000e+01, 2.0160000000e+01, 5.7200000000e+00,
                2.1820000000e+00, 4.3270000000e-01, 1.5700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5010000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0700000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.4754100000e-03, 3.4170800000e-02, 1.4425000000e-01,
                           3.5392800000e-01, 4.5908500000e-01, 2.0638300000e-01,
                           1.0214100000e-02, -6.0312200000e-05},
                 doubles_t{3.9970000000e+02, 9.4190000000e+01, 2.9750000000e+01,
                           1.0770000000e+01, 4.1190000000e+00, 1.6250000000e+00,
                           4.7260000000e-01, 1.4070000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.1625100000e-03, -8.6566400000e-03, -3.9088600000e-02,
                -9.3462500000e-02, -1.4799400000e-01, 3.0190400000e-02,
                5.6157300000e-01, 5.3477600000e-01},
      doubles_t{3.9970000000e+02, 9.4190000000e+01, 2.9750000000e+01,
                1.0770000000e+01, 4.1190000000e+00, 1.6250000000e+00,
                4.7260000000e-01, 1.4070000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4070000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5030000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.9900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2110000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5200000000e-01}));
    return abs_t(name, 16, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pwcvdz_16

} // namespace chemcache
