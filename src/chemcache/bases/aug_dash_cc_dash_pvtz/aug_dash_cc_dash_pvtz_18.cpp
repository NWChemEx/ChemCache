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

#include "aug_dash_cc_dash_pvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvtz_18() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.5582800000e-05, 3.5410800000e-04, 1.8579700000e-03,
                7.7685100000e-03, 2.7423200000e-02, 8.2383600000e-02,
                2.0123000000e-01, 3.5678100000e-01, 3.4956300000e-01,
                1.1826600000e-01, 5.6019000000e-03, 4.8347300000e-04,
                -1.4808500000e-04, 2.9202500000e-05, -2.3160400000e-05},
      doubles_t{5.4500000000e+05, 8.1640000000e+04, 1.8580000000e+04,
                5.2610000000e+03, 1.7170000000e+03, 6.1990000000e+02,
                2.4160000000e+02, 9.9790000000e+01, 4.3150000000e+01,
                1.9140000000e+01, 7.4880000000e+00, 3.2050000000e+00,
                1.1960000000e+00, 5.2040000000e-01, 1.9540000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.2955100000e-05, -1.0042800000e-04, -5.2958300000e-04,
                -2.2139600000e-03, -7.9684500000e-03, -2.4580300000e-02,
                -6.5779800000e-02, -1.3794200000e-01, -2.0163000000e-01,
                -4.1283400000e-02, 4.8468000000e-01, 5.7922400000e-01,
                8.7908300000e-02, -7.2755300000e-03, 2.3288400000e-03},
      doubles_t{5.4500000000e+05, 8.1640000000e+04, 1.8580000000e+04,
                5.2610000000e+03, 1.7170000000e+03, 6.1990000000e+02,
                2.4160000000e+02, 9.9790000000e+01, 4.3150000000e+01,
                1.9140000000e+01, 7.4880000000e+00, 3.2050000000e+00,
                1.1960000000e+00, 5.2040000000e-01, 1.9540000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1960000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.0499000000e-06, 3.1369100000e-05, 1.6564600000e-04,
                6.9166200000e-04, 2.4979000000e-03, 7.7107400000e-03,
                2.0871400000e-02, 4.4396500000e-02, 6.8022400000e-02,
                1.4135000000e-02, -2.0748900000e-01, -4.2504500000e-01,
                7.0485500000e-02, 7.3362700000e-01, 3.9600500000e-01},
      doubles_t{5.4500000000e+05, 8.1640000000e+04, 1.8580000000e+04,
                5.2610000000e+03, 1.7170000000e+03, 6.1990000000e+02,
                2.4160000000e+02, 9.9790000000e+01, 4.3150000000e+01,
                1.9140000000e+01, 7.4880000000e+00, 3.2050000000e+00,
                1.1960000000e+00, 5.2040000000e-01, 1.9540000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9540000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8500000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.3697600000e-03, 1.9019900000e-02, 8.8080700000e-02,
                2.5637700000e-01, 4.3871100000e-01, 3.4756900000e-01,
                5.6674000000e-02, -5.2388200000e-03, 1.6437600000e-03},
      doubles_t{7.6180000000e+02, 1.8020000000e+02, 5.7500000000e+01,
                2.1240000000e+01, 8.3880000000e+00, 3.4160000000e+00,
                1.2060000000e+00, 4.5230000000e-01, 1.5450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2060000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-6.6721100000e-04, -5.3271700000e-03, -2.5549400000e-02,
                -7.5719700000e-02, -1.4113300000e-01, -9.3276800000e-02,
                2.8287200000e-01, 5.6245000000e-01, 3.2505900000e-01},
      doubles_t{7.6180000000e+02, 1.8020000000e+02, 5.7500000000e+01,
                2.1240000000e+01, 8.3880000000e+00, 3.4160000000e+00,
                1.2060000000e+00, 4.5230000000e-01, 1.5450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5450000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2540000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.9000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0600000000e-01}));
    return abs_t(name, 18, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvtz_18

} // namespace chemcache
