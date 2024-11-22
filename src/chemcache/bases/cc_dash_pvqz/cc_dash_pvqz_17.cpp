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

#include "cc_dash_pvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvqz_17() {
    // Basis Set name and origin point
    std::string name("cc-pvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.3168800000e-05, 1.8015400000e-04, 9.4778200000e-04,
                4.0013900000e-03, 1.4462900000e-02, 4.5658600000e-02,
                1.2324800000e-01, 2.6436900000e-01, 3.8298900000e-01,
                2.7093400000e-01, 4.7140400000e-02, -3.7176600000e-03,
                2.1915800000e-03, -8.4985000000e-04, 4.2466400000e-04,
                -1.2001800000e-04},
      doubles_t{8.3490000000e+05, 1.2500000000e+05, 2.8430000000e+04,
                8.0330000000e+03, 2.6080000000e+03, 9.3390000000e+02,
                3.6000000000e+02, 1.4700000000e+02, 6.2880000000e+01,
                2.7600000000e+01, 1.1080000000e+01, 5.0750000000e+00,
                2.2780000000e+00, 7.7750000000e-01, 3.5270000000e-01,
                1.4310000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.4964900000e-06, -5.0489500000e-05, -2.6611300000e-04,
                -1.1249900000e-03, -4.1049700000e-03, -1.3198700000e-02,
                -3.7534200000e-02, -8.9723300000e-02, -1.6767100000e-01,
                -1.7476300000e-01, 1.1490900000e-01, 5.6361800000e-01,
                4.4160600000e-01, 4.2669700000e-02, -6.6722900000e-03,
                1.9072900000e-03},
      doubles_t{8.3490000000e+05, 1.2500000000e+05, 2.8430000000e+04,
                8.0330000000e+03, 2.6080000000e+03, 9.3390000000e+02,
                3.6000000000e+02, 1.4700000000e+02, 6.2880000000e+01,
                2.7600000000e+01, 1.1080000000e+01, 5.0750000000e+00,
                2.2780000000e+00, 7.7750000000e-01, 3.5270000000e-01,
                1.4310000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7750000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5270000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9664500000e-06, 1.5262000000e-05, 8.0608600000e-05,
                3.3996000000e-04, 1.2455100000e-03, 3.9961200000e-03,
                1.1475100000e-02, 2.7550400000e-02, 5.3291700000e-02,
                5.7124600000e-02, -3.9520100000e-02, -2.6434300000e-01,
                -3.4929100000e-01, 2.6967100000e-01, 6.7607300000e-01,
                2.8467900000e-01},
      doubles_t{8.3490000000e+05, 1.2500000000e+05, 2.8430000000e+04,
                8.0330000000e+03, 2.6080000000e+03, 9.3390000000e+02,
                3.6000000000e+02, 1.4700000000e+02, 6.2880000000e+01,
                2.7600000000e+01, 1.1080000000e+01, 5.0750000000e+00,
                2.2780000000e+00, 7.7750000000e-01, 3.5270000000e-01,
                1.4310000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4310000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.7403900000e-04, 4.0641200000e-03, 2.1335500000e-02,
                           7.9461100000e-02, 2.0892700000e-01, 3.6494500000e-01,
                           3.7172500000e-01, 1.4629200000e-01, 1.0790600000e-02,
                           1.1700400000e-03, 3.3940800000e-04},
                 doubles_t{1.7030000000e+03, 4.0360000000e+02, 1.3030000000e+02,
                           4.9050000000e+01, 2.0260000000e+01, 8.7870000000e+00,
                           3.9190000000e+00, 1.7650000000e+00, 7.2070000000e-01,
                           2.8390000000e-01, 1.0600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.2070000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.2826600000e-04, -1.0935600000e-03, -5.8342900000e-03,
                -2.1925800000e-02, -6.0138500000e-02, -1.0692900000e-01,
                -1.2245400000e-01, 3.8361900000e-02, 3.8525600000e-01,
                5.0726500000e-01, 2.2721800000e-01},
      doubles_t{1.7030000000e+03, 4.0360000000e+02, 1.3030000000e+02,
                4.9050000000e+01, 2.0260000000e+01, 8.7870000000e+00,
                3.9190000000e+00, 1.7650000000e+00, 7.2070000000e-01,
                2.8390000000e-01, 1.0600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8390000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5510000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0890000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{8.2700000000e-01}));
    return abs_t(name, 17, r0, shells.begin(), shells.end());
} // cc_dash_pvqz_17

} // namespace chemcache
