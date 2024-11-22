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

#include "cc_dash_pv_oparen_q_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_q_plus_d_cparen_z_12() {
    // Basis Set name and origin point
    std::string name("cc-pv(q+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{3.0960800000e-05, 2.4095400000e-04, 1.2666000000e-03,
                5.3335900000e-03, 1.9077000000e-02, 5.8805800000e-02,
                1.5145400000e-01, 3.0071600000e-01, 3.8114900000e-01,
                2.1358400000e-01, 2.3121000000e-02, -2.3075700000e-03,
                1.2890000000e-03, -4.2638700000e-04, 3.5432000000e-04,
                -1.1258800000e-04},
      doubles_t{3.2760000000e+05, 4.9050000000e+04, 1.1150000000e+04,
                3.1520000000e+03, 1.0250000000e+03, 3.6880000000e+02,
                1.4320000000e+02, 5.8960000000e+01, 2.5400000000e+01,
                1.1150000000e+01, 4.0040000000e+00, 1.7010000000e+00,
                7.0600000000e-01, 1.4100000000e-01, 6.8080000000e-02,
                3.0630000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.8317300000e-06, -6.0793500000e-05, -3.2119700000e-04,
                -1.3495500000e-03, -4.9057000000e-03, -1.5356100000e-02,
                -4.2340900000e-02, -9.4060300000e-02, -1.6342500000e-01,
                -1.2475400000e-01, 2.3562300000e-01, 5.7756300000e-01,
                3.3523200000e-01, 1.5633800000e-02, -8.3779900000e-03,
                2.2305600000e-03},
      doubles_t{3.2760000000e+05, 4.9050000000e+04, 1.1150000000e+04,
                3.1520000000e+03, 1.0250000000e+03, 3.6880000000e+02,
                1.4320000000e+02, 5.8960000000e+01, 2.5400000000e+01,
                1.1150000000e+01, 4.0040000000e+00, 1.7010000000e+00,
                7.0600000000e-01, 1.4100000000e-01, 6.8080000000e-02,
                3.0630000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8080000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.5090800000e-06, 1.1713400000e-05, 6.1898000000e-05,
                2.6008800000e-04, 9.4621800000e-04, 2.9659500000e-03,
                8.2124500000e-03, 1.8397700000e-02, 3.2665700000e-02,
                2.5731500000e-02, -5.3535100000e-02, -1.5689500000e-01,
                -2.0665900000e-01, 3.2429800000e-01, 5.5261100000e-01,
                2.6035200000e-01},
      doubles_t{3.2760000000e+05, 4.9050000000e+04, 1.1150000000e+04,
                3.1520000000e+03, 1.0250000000e+03, 3.6880000000e+02,
                1.4320000000e+02, 5.8960000000e+01, 2.5400000000e+01,
                1.1150000000e+01, 4.0040000000e+00, 1.7010000000e+00,
                7.0600000000e-01, 1.4100000000e-01, 6.8080000000e-02,
                3.0630000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0630000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{8.3396900000e-04, 6.8921500000e-03, 3.3787400000e-02,
                1.1440100000e-01, 2.5951400000e-01, 3.8509500000e-01,
                3.3537300000e-01, 1.1064100000e-01, -1.2131500000e-02,
                7.0056200000e-03, -1.5133500000e-03, 4.4892500000e-04},
      doubles_t{5.3960000000e+02, 1.2790000000e+02, 4.1020000000e+01,
                1.5250000000e+01, 6.1660000000e+00, 2.5610000000e+00,
                1.0600000000e+00, 4.1760000000e-01, 2.6900000000e-01,
                1.2230000000e-01, 5.4760000000e-02, 2.3880000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2230000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.3207600000e-04, -1.0953800000e-03, -5.3949500000e-03,
                -1.8557200000e-02, -4.2737500000e-02, -6.4768400000e-02,
                -6.2781800000e-02, -2.4491200000e-02, 1.0476100000e-01,
                3.5890900000e-01, 4.9250100000e-01, 1.8528000000e-01},
      doubles_t{5.3960000000e+02, 1.2790000000e+02, 4.1020000000e+01,
                1.5250000000e+01, 6.1660000000e+00, 2.5610000000e+00,
                1.0600000000e+00, 4.1760000000e-01, 2.6900000000e-01,
                1.2230000000e-01, 5.4760000000e-02, 2.3880000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4760000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3880000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9950000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7290000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0280000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6220000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7910000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0630000000e-01}));
    return abs_t(name, 12, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_q_plus_d_cparen_z_12

} // namespace chemcache
