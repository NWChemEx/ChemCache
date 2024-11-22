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

#include "aug_dash_cc_dash_pwcv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pwcv5z_14() {
    // Basis Set name and origin point
    std::string name("aug-cc-pwcv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0371200000e-06,  1.5839400000e-05,  8.3359000000e-05,
                3.5136100000e-04,  1.2766000000e-03,  4.1519100000e-03,
                1.2303000000e-02,  3.3310200000e-02,  8.0984500000e-02,
                1.7029000000e-01,  2.8687900000e-01,  3.3034000000e-01,
                1.9660200000e-01,  3.5453500000e-02,  -5.3520400000e-04,
                1.6146500000e-03,  -3.7273900000e-04, 1.4622600000e-04,
                -7.8944600000e-05, 1.9280100000e-05},
      doubles_t{3.9480000000e+06, 5.9110000000e+05, 1.3450000000e+05,
                3.8120000000e+04, 1.2460000000e+04, 4.5040000000e+03,
                1.7580000000e+03, 7.2910000000e+02, 3.1800000000e+02,
                1.4460000000e+02, 6.7970000000e+01, 3.2820000000e+01,
                1.6030000000e+01, 7.3960000000e+00, 3.6610000000e+00,
                1.8230000000e+00, 9.1470000000e-01, 3.3930000000e-01,
                1.5000000000e-01, 6.4380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2820000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6030000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3960000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6610000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.4208500000e-07, -4.2167700000e-06, -2.2181300000e-05,
                -9.3602800000e-05, -3.4011600000e-04, -1.1106100000e-03,
                -3.3087800000e-03, -9.1160200000e-03, -2.2879000000e-02,
                -5.1711900000e-02, -9.9909100000e-02, -1.5274700000e-01,
                -1.2750800000e-01, 9.4696300000e-02,  4.1403600000e-01,
                4.6793400000e-01,  1.7392700000e-01,  8.4389500000e-03,
                -9.9806700000e-04, 3.6209600000e-04},
      doubles_t{3.9480000000e+06, 5.9110000000e+05, 1.3450000000e+05,
                3.8120000000e+04, 1.2460000000e+04, 4.5040000000e+03,
                1.7580000000e+03, 7.2910000000e+02, 3.1800000000e+02,
                1.4460000000e+02, 6.7970000000e+01, 3.2820000000e+01,
                1.6030000000e+01, 7.3960000000e+00, 3.6610000000e+00,
                1.8230000000e+00, 9.1470000000e-01, 3.3930000000e-01,
                1.5000000000e-01, 6.4380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8230000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1470000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3930000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1810000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1600000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5630000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3500000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.0120600000e-04, 1.7493700000e-03, 9.4814100000e-03,
                3.7231300000e-02, 1.1076300000e-01, 2.3793300000e-01,
                3.5369100000e-01, 3.2883900000e-01, 1.3237300000e-01,
                1.0330000000e-02, -1.5030500000e-04, 2.6580900000e-04},
      doubles_t{1.7800000000e+03, 4.2180000000e+02, 1.3670000000e+02,
                5.1810000000e+01, 2.1600000000e+01, 9.5630000000e+00,
                4.3500000000e+00, 2.0060000000e+00, 9.2050000000e-01,
                3.5000000000e-01, 1.3810000000e-01, 5.3380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0060000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.2050000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3810000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3380000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9732000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3420000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0570000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7340000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0820000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.1700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3381000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1010000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9450000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3500000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1435000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1610000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2250000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2300000000e-01}));
    return abs_t(name, 14, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pwcv5z_14

} // namespace chemcache
