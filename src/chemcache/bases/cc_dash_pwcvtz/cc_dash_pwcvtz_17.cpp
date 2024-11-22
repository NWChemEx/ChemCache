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

#include "cc_dash_pwcvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pwcvtz_17() {
    // Basis Set name and origin point
    std::string name("cc-pwcvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.9297000000e-05, 3.8302900000e-04, 2.0085400000e-03,
                8.3855800000e-03, 2.9470300000e-02, 8.7832500000e-02,
                2.1147300000e-01, 3.6536400000e-01, 3.4088400000e-01,
                1.0213300000e-01, 3.1167500000e-03, 1.0575100000e-03,
                -3.7800000000e-04, 1.5613600000e-04, -5.1412600000e-05},
      doubles_t{4.5610000000e+05, 6.8330000000e+04, 1.5550000000e+04,
                4.4050000000e+03, 1.4390000000e+03, 5.2040000000e+02,
                2.0310000000e+02, 8.3960000000e+01, 3.6200000000e+01,
                1.5830000000e+01, 6.3340000000e+00, 2.6940000000e+00,
                9.7680000000e-01, 4.3130000000e-01, 1.6250000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.3830400000e-05, -1.0727900000e-04, -5.6508300000e-04,
                -2.3613500000e-03, -8.4588600000e-03, -2.5963800000e-02,
                -6.8636200000e-02, -1.4187400000e-01, -1.9931900000e-01,
                -1.9566200000e-02, 4.9974100000e-01, 5.6373600000e-01,
                7.9032500000e-02, -8.3509100000e-03, 2.3245600000e-03},
      doubles_t{4.5610000000e+05, 6.8330000000e+04, 1.5550000000e+04,
                4.4050000000e+03, 1.4390000000e+03, 5.2040000000e+02,
                2.0310000000e+02, 8.3960000000e+01, 3.6200000000e+01,
                1.5830000000e+01, 6.3340000000e+00, 2.6940000000e+00,
                9.7680000000e-01, 4.3130000000e-01, 1.6250000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.7680000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.1854600000e-06, 3.2439500000e-05, 1.7110500000e-04,
                7.1417600000e-04, 2.5670500000e-03, 7.8855200000e-03,
                2.1086700000e-02, 4.4226400000e-02, 6.5167000000e-02,
                6.0301200000e-03, -2.0649500000e-01, -4.0587100000e-01,
                7.5955800000e-02, 7.2566100000e-01, 3.9442300000e-01},
      doubles_t{4.5610000000e+05, 6.8330000000e+04, 1.5550000000e+04,
                4.4050000000e+03, 1.4390000000e+03, 5.2040000000e+02,
                2.0310000000e+02, 8.3960000000e+01, 3.6200000000e+01,
                1.5830000000e+01, 6.3340000000e+00, 2.6940000000e+00,
                9.7680000000e-01, 4.3130000000e-01, 1.6250000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6250000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3940000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5620000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6930000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.4044800000e-03, 1.9214800000e-02, 8.8509700000e-02,
                2.5602000000e-01, 4.3692700000e-01, 3.5033400000e-01,
                5.8549500000e-02, -4.5842300000e-03, 2.2697000000e-03},
      doubles_t{6.6330000000e+02, 1.5680000000e+02, 4.9980000000e+01,
                1.8420000000e+01, 7.2400000000e+00, 2.9220000000e+00,
                1.0220000000e+00, 3.8180000000e-01, 1.3010000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0220000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-6.5214500000e-04, -5.1944500000e-03, -2.4693800000e-02,
                -7.2816700000e-02, -1.3403000000e-01, -9.4774200000e-02,
                2.6228900000e-01, 5.6466700000e-01, 3.4125000000e-01},
      doubles_t{6.6330000000e+02, 1.5680000000e+02, 4.9980000000e+01,
                1.8420000000e+01, 7.2400000000e+00, 2.9220000000e+00,
                1.0220000000e+00, 3.8180000000e-01, 1.3010000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3010000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1170000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3506000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5370000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0460000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0140000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0600000000e-01}));
    return abs_t(name, 17, r0, shells.begin(), shells.end());
} // cc_dash_pwcvtz_17

} // namespace chemcache
