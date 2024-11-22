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

#include "aug_dash_cc_dash_pvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvqz_6() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.1000000000e-05, 7.0400000000e-04, 3.6930000000e-03,
                1.5360000000e-02, 5.2929000000e-02, 1.4704300000e-01,
                3.0563100000e-01, 3.9934500000e-01, 2.1705100000e-01,
                1.5894000000e-02, -3.0840000000e-03, 9.7800000000e-04},
      doubles_t{3.3980000000e+04, 5.0890000000e+03, 1.1570000000e+03,
                3.2660000000e+02, 1.0610000000e+02, 3.8110000000e+01,
                1.4750000000e+01, 6.0350000000e+00, 2.5300000000e+00,
                7.3550000000e-01, 2.9050000000e-01, 1.1110000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3550000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9050000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.9000000000e-05, -1.5100000000e-04, -7.8500000000e-04,
                -3.3240000000e-03, -1.1512000000e-02, -3.4160000000e-02,
                -7.7173000000e-02, -1.4149300000e-01, -1.1801900000e-01,
                2.7380600000e-01, 5.8651000000e-01, 2.8543000000e-01},
      doubles_t{3.3980000000e+04, 5.0890000000e+03, 1.1570000000e+03,
                3.2660000000e+02, 1.0610000000e+02, 3.8110000000e+01,
                1.4750000000e+01, 6.0350000000e+00, 2.5300000000e+00,
                7.3550000000e-01, 2.9050000000e-01, 1.1110000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1110000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1450000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.1320000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.3780000000e-03, 3.6132000000e-02, 1.4249300000e-01,
                3.4215000000e-01, 4.6386400000e-01, 2.5002800000e-01},
      doubles_t{3.4510000000e+01, 7.9150000000e+00, 2.3680000000e+00,
                8.1320000000e-01, 2.8900000000e-01, 1.0070000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0070000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2180000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8480000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4190000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0110000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2400000000e-01}));
    return abs_t(name, 6, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvqz_6

} // namespace chemcache
