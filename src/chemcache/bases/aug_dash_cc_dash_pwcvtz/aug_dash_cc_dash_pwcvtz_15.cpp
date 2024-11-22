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

#include "aug_dash_cc_dash_pwcvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pwcvtz_15() {
    // Basis Set name and origin point
    std::string name("aug-cc-pwcvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.7696000000e-05, 4.4829600000e-04, 2.3493900000e-03,
                9.7826500000e-03, 3.4146700000e-02, 1.0020400000e-01,
                2.3437200000e-01, 3.8243400000e-01, 3.1808800000e-01,
                7.0778800000e-02, -1.8179900000e-03, 2.1618000000e-03,
                -8.3474200000e-04, 4.3229700000e-04, -1.1425100000e-04},
      doubles_t{3.1240000000e+05, 4.6800000000e+04, 1.0650000000e+04,
                3.0180000000e+03, 9.8680000000e+02, 3.5740000000e+02,
                1.3960000000e+02, 5.7630000000e+01, 2.4600000000e+01,
                1.0120000000e+01, 4.2830000000e+00, 1.8050000000e+00,
                6.1580000000e-01, 2.7820000000e-01, 1.0550000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.5670900000e-05, -1.2172400000e-04, -6.3967200000e-04,
                -2.6742600000e-03, -9.4983100000e-03, -2.8934900000e-02,
                -7.4512100000e-02, -1.4993800000e-01, -1.8946700000e-01,
                3.6327000000e-02, 5.2881600000e-01, 5.1911500000e-01,
                6.0554700000e-02, -9.2569500000e-03, 2.1037200000e-03},
      doubles_t{3.1240000000e+05, 4.6800000000e+04, 1.0650000000e+04,
                3.0180000000e+03, 9.8680000000e+02, 3.5740000000e+02,
                1.3960000000e+02, 5.7630000000e+01, 2.4600000000e+01,
                1.0120000000e+01, 4.2830000000e+00, 1.8050000000e+00,
                6.1580000000e-01, 2.7820000000e-01, 1.0550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1580000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.3063100000e-06, 3.3419400000e-05, 1.7588500000e-04,
                7.3434000000e-04, 2.6177500000e-03, 7.9785200000e-03,
                2.0794000000e-02, 4.2444600000e-02, 5.6343600000e-02,
                -1.2735800000e-02, -1.9649500000e-01, -3.5355500000e-01,
                7.4140700000e-02, 7.0091200000e-01, 4.0473900000e-01},
      doubles_t{3.1240000000e+05, 4.6800000000e+04, 1.0650000000e+04,
                3.0180000000e+03, 9.8680000000e+02, 3.5740000000e+02,
                1.3960000000e+02, 5.7630000000e+01, 2.4600000000e+01,
                1.0120000000e+01, 4.2830000000e+00, 1.8050000000e+00,
                6.1580000000e-01, 2.7820000000e-01, 1.0550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0069000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7300000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0900000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.3372800000e-03, 1.8541000000e-02, 8.4969300000e-02,
                2.4461500000e-01, 4.2276600000e-01, 3.6843900000e-01,
                7.7273400000e-02, -3.7900500000e-03, 1.5993900000e-03},
      doubles_t{5.0490000000e+02, 1.1940000000e+02, 3.7960000000e+01,
                1.3950000000e+01, 5.4570000000e+00, 2.1770000000e+00,
                8.0100000000e-01, 2.8770000000e-01, 9.7140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0100000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-5.5523600000e-04, -4.4591300000e-03, -2.0635000000e-02,
                -6.1769400000e-02, -1.0892400000e-01, -1.0559900000e-01,
                1.5348300000e-01, 5.7698100000e-01, 4.2243900000e-01},
      doubles_t{5.0490000000e+02, 1.1940000000e+02, 3.7960000000e+01,
                1.3950000000e+01, 5.4570000000e+00, 2.1770000000e+00,
                8.0100000000e-01, 2.8770000000e-01, 9.7140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.7140000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4010000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.0970000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3280000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1220000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6500000000e-01}));
    return abs_t(name, 15, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pwcvtz_15

} // namespace chemcache
