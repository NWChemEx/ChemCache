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

abs_t aug_dash_cc_dash_pvtz_13() {
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
      doubles_t{6.7883600000e-05, 5.2714900000e-04, 2.7620300000e-03,
                1.1472800000e-02, 3.9818800000e-02, 1.1504000000e-01,
                2.6088700000e-01, 3.9638600000e-01, 2.8459700000e-01,
                4.4458300000e-02, -4.8983800000e-03, 2.6125300000e-03,
                -1.0891500000e-03, 7.2206800000e-04, -1.7848500000e-04},
      doubles_t{2.0550000000e+05, 3.0780000000e+04, 7.0060000000e+03,
                1.9850000000e+03, 6.4910000000e+02, 2.3500000000e+02,
                9.1620000000e+01, 3.7670000000e+01, 1.5910000000e+01,
                5.8500000000e+00, 2.5420000000e+00, 1.0570000000e+00,
                2.9310000000e-01, 1.4550000000e-01, 5.6500000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7637700000e-05, -1.3719500000e-04, -7.1891000000e-04,
                -3.0114600000e-03, -1.0601400000e-02, -3.2134500000e-02,
                -8.0315600000e-02, -1.5679400000e-01, -1.6837600000e-01,
                1.2687900000e-01, 5.6149400000e-01, 4.3661300000e-01,
                3.5729300000e-02, -1.1456300000e-02, 2.2020100000e-03},
      doubles_t{2.0550000000e+05, 3.0780000000e+04, 7.0060000000e+03,
                1.9850000000e+03, 6.4910000000e+02, 2.3500000000e+02,
                9.1620000000e+01, 3.7670000000e+01, 1.5910000000e+01,
                5.8500000000e+00, 2.5420000000e+00, 1.0570000000e+00,
                2.9310000000e-01, 1.4550000000e-01, 5.6500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9310000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.0731500000e-06, 3.1656600000e-05, 1.6611600000e-04,
                6.9499200000e-04, 2.4551100000e-03, 7.4459800000e-03,
                1.8825300000e-02, 3.7277200000e-02, 4.1949600000e-02,
                -3.5437500000e-02, -1.7513200000e-01, -2.7620300000e-01,
                1.0872900000e-01, 6.5280900000e-01, 3.9458700000e-01},
      doubles_t{2.0550000000e+05, 3.0780000000e+04, 7.0060000000e+03,
                1.9850000000e+03, 6.4910000000e+02, 2.3500000000e+02,
                9.1620000000e+01, 3.7670000000e+01, 1.5910000000e+01,
                5.8500000000e+00, 2.5420000000e+00, 1.0570000000e+00,
                2.9310000000e-01, 1.4550000000e-01, 5.6500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.6500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2100000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.6278600000e-03, 1.3068700000e-02, 6.1234100000e-02,
                1.8787000000e-01, 3.6045200000e-01, 4.0845400000e-01,
                1.8864000000e-01, 9.7651400000e-03, -1.1505700000e-03},
      doubles_t{4.4440000000e+02, 1.0510000000e+02, 3.3470000000e+01,
                1.2330000000e+01, 4.8690000000e+00, 1.9610000000e+00,
                7.8340000000e-01, 1.8880000000e-01, 5.5570000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8340000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-2.8634100000e-04, -2.4230800000e-03, -1.0865800000e-02,
                -3.6430700000e-02, -6.4107400000e-02, -9.7223900000e-02,
                1.4743700000e-02, 5.0344800000e-01, 5.9798400000e-01},
      doubles_t{4.4440000000e+02, 1.0510000000e+02, 3.3470000000e+01,
                1.2330000000e+01, 4.8690000000e+00, 1.9610000000e+00,
                7.8340000000e-01, 1.8880000000e-01, 5.5570000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5570000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5800000000e-02}));
    return abs_t(name, 13, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvtz_13

} // namespace chemcache
