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

#include "aug_dash_cc_dash_pv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv5z_7() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{2.5000000000e-05, 1.9700000000e-04, 1.0320000000e-03,
                           4.3250000000e-03, 1.5380000000e-02, 4.6867000000e-02,
                           1.2011600000e-01, 2.4569500000e-01, 3.6137900000e-01,
                           2.8728300000e-01, 7.0171000000e-02, 1.8310000000e-03,
                           8.3500000000e-04, -6.0000000000e-06},
                 doubles_t{1.2920000000e+05, 1.9350000000e+04, 4.4040000000e+03,
                           1.2480000000e+03, 4.0800000000e+02, 1.4820000000e+02,
                           5.8500000000e+01, 2.4590000000e+01, 1.0810000000e+01,
                           4.8820000000e+00, 2.1950000000e+00, 8.7150000000e-01,
                           3.5040000000e-01, 1.3970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1950000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7150000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5040000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.0000000000e-06, -4.3000000000e-05, -2.2700000000e-04,
                -9.5800000000e-04, -3.4160000000e-03, -1.0667000000e-02,
                -2.8279000000e-02, -6.4020000000e-02, -1.1393200000e-01,
                -1.4699500000e-01, -7.2510000000e-03, 3.6618300000e-01,
                5.4790800000e-01, 2.1664500000e-01},
      doubles_t{1.2920000000e+05, 1.9350000000e+04, 4.4040000000e+03,
                1.2480000000e+03, 4.0800000000e+02, 1.4820000000e+02,
                5.8500000000e+01, 2.4590000000e+01, 1.0810000000e+01,
                4.8820000000e+00, 2.1950000000e+00, 8.7150000000e-01,
                3.5040000000e-01, 1.3970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3970000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1800000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5870000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5330000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{8.9200000000e-04, 7.0820000000e-03, 3.2816000000e-02,
                           1.0820900000e-01, 2.4809400000e-01, 3.7451300000e-01,
                           3.4841400000e-01, 1.2834000000e-01},
                 doubles_t{1.4700000000e+02, 3.4760000000e+01, 1.1000000000e+01,
                           3.9950000000e+00, 1.5870000000e+00, 6.5330000000e-01,
                           2.6860000000e-01, 1.0670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6860000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6900000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6470000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8130000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.7100000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9420000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2040000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5110000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{9.4200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7680000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8800000000e-01}));
    return abs_t(name, 7, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv5z_7

} // namespace chemcache
