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

abs_t cc_dash_pvdz_31() {
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
      doubles_t{2.0680000000e-04, 1.6047000000e-03, 8.3402000000e-03,
                3.4024800000e-02, 1.1116990000e-01, 2.7539300000e-01,
                4.2126280000e-01, 2.7389060000e-01, 2.8372000000e-02,
                -6.2931000000e-03, 2.0606000000e-03, -9.2690000000e-04,
                2.2730000000e-04, -1.0630000000e-04},
      doubles_t{4.8513000000e+05, 7.2719000000e+04, 1.6552000000e+04,
                4.6878000000e+03, 1.5291000000e+03, 5.5181000000e+02,
                2.1518000000e+02, 8.8174000000e+01, 2.7154000000e+01,
                1.1503000000e+01, 3.3018000000e+00, 1.3314000000e+00,
                1.9316000000e-01, 7.0895000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.4300000000e-05, -4.9540000000e-04, -2.6208000000e-03,
                -1.0683900000e-02, -3.7412300000e-02, -1.0096360000e-01,
                -2.1451410000e-01, -1.7522970000e-01, 4.8315990000e-01,
                6.3236770000e-01, 6.8494200000e-02, -1.1871200000e-02,
                2.6652000000e-03, -1.2251000000e-03},
      doubles_t{4.8513000000e+05, 7.2719000000e+04, 1.6552000000e+04,
                4.6878000000e+03, 1.5291000000e+03, 5.5181000000e+02,
                2.1518000000e+02, 8.8174000000e+01, 2.7154000000e+01,
                1.1503000000e+01, 3.3018000000e+00, 1.3314000000e+00,
                1.9316000000e-01, 7.0895000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4500000000e-05, 1.8950000000e-04, 9.9640000000e-04,
                4.1082000000e-03, 1.4293800000e-02, 3.9803400000e-02,
                8.5594000000e-02, 7.9630500000e-02, -2.9391070000e-01,
                -5.2639140000e-01, 5.8642490000e-01, 6.7263470000e-01,
                2.7612300000e-02, -9.3651000000e-03},
      doubles_t{4.8513000000e+05, 7.2719000000e+04, 1.6552000000e+04,
                4.6878000000e+03, 1.5291000000e+03, 5.5181000000e+02,
                2.1518000000e+02, 8.8174000000e+01, 2.7154000000e+01,
                1.1503000000e+01, 3.3018000000e+00, 1.3314000000e+00,
                1.9316000000e-01, 7.0895000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.7000000000e-06, -4.4000000000e-05, -2.3050000000e-04,
                -9.5440000000e-04, -3.3055000000e-03, -9.2888000000e-03,
                -1.9864400000e-02, -1.9088800000e-02, 7.3235600000e-02,
                1.3415260000e-01, -1.8319290000e-01, -3.5713080000e-01,
                6.2460130000e-01, 5.2384300000e-01},
      doubles_t{4.8513000000e+05, 7.2719000000e+04, 1.6552000000e+04,
                4.6878000000e+03, 1.5291000000e+03, 5.5181000000e+02,
                2.1518000000e+02, 8.8174000000e+01, 2.7154000000e+01,
                1.1503000000e+01, 3.3018000000e+00, 1.3314000000e+00,
                1.9316000000e-01, 7.0895000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0895000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.5260000000e-03, 1.2748600000e-02, 6.3374200000e-02,
                2.0657750000e-01, 4.0929630000e-01, 3.9191830000e-01,
                1.0294410000e-01, -7.2030000000e-04, 2.0950000000e-03,
                -3.2900000000e-04, 1.1620000000e-04},
      doubles_t{3.2486000000e+03, 7.6997000000e+02, 2.4820000000e+02,
                9.3364000000e+01, 3.8251000000e+01, 1.6422000000e+01,
                6.7918000000e+00, 2.8336000000e+00, 1.1062000000e+00,
                2.2250000000e-01, 6.1772000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-5.8030000000e-04, -4.8647000000e-03, -2.4839400000e-02,
                -8.4175900000e-02, -1.8008850000e-01, -1.5855550000e-01,
                2.3553760000e-01, 5.8205870000e-01, 3.3666190000e-01,
                1.7191200000e-02, -3.3265000000e-03},
      doubles_t{3.2486000000e+03, 7.6997000000e+02, 2.4820000000e+02,
                9.3364000000e+01, 3.8251000000e+01, 1.6422000000e+01,
                6.7918000000e+00, 2.8336000000e+00, 1.1062000000e+00,
                2.2250000000e-01, 6.1772000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{9.5000000000e-05, 7.8320000000e-04, 4.0855000000e-03,
                1.3598700000e-02, 3.0269500000e-02, 2.4179000000e-02,
                -4.2377700000e-02, -1.2656610000e-01, -4.9944400000e-02,
                4.4941990000e-01, 6.7188990000e-01},
      doubles_t{3.2486000000e+03, 7.6997000000e+02, 2.4820000000e+02,
                9.3364000000e+01, 3.8251000000e+01, 1.6422000000e+01,
                6.7918000000e+00, 2.8336000000e+00, 1.1062000000e+00,
                2.2250000000e-01, 6.1772000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1772000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.7382500000e-02, 1.5108050000e-01, 3.7492170000e-01,
                           4.7507990000e-01, 2.9827500000e-01},
                 doubles_t{6.5337000000e+01, 1.8497000000e+01, 6.3150000000e+00,
                           2.1635000000e+00, 6.6675000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8840000000e-01}));
    return abs_t(name, 31, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_31

} // namespace chemcache
