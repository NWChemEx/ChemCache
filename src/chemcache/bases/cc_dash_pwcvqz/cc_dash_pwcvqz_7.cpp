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

#include "cc_dash_pwcvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pwcvqz_7() {
    // Basis Set name and origin point
    std::string name("cc-pwcvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4917000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0735000000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.2000000000e-05, 7.1700000000e-04, 3.7490000000e-03,
                1.5532000000e-02, 5.3146000000e-02, 1.4678700000e-01,
                3.0466300000e-01, 3.9768400000e-01, 2.1764100000e-01,
                1.6963000000e-02, -2.7450000000e-03, 9.5300000000e-04},
      doubles_t{4.5840000000e+04, 6.8680000000e+03, 1.5630000000e+03,
                4.4240000000e+02, 1.4430000000e+02, 5.2180000000e+01,
                2.0340000000e+01, 8.3810000000e+00, 3.5290000000e+00,
                1.0540000000e+00, 4.1180000000e-01, 1.5520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0540000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1180000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.0000000000e-05, -1.5900000000e-04, -8.2400000000e-04,
                -3.4780000000e-03, -1.1966000000e-02, -3.5388000000e-02,
                -8.0077000000e-02, -1.4672200000e-01, -1.1636000000e-01,
                2.7991900000e-01, 5.8548100000e-01, 2.8402800000e-01},
      doubles_t{4.5840000000e+04, 6.8680000000e+03, 1.5630000000e+03,
                4.4240000000e+02, 1.4430000000e+02, 5.2180000000e+01,
                2.0340000000e+01, 8.3810000000e+00, 3.5290000000e+00,
                1.0540000000e+00, 4.1180000000e-01, 1.5520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5520000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8290000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.4185000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5190000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5540000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1820000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.5330000000e-03, 3.7962000000e-02, 1.4902800000e-01,
                3.4892200000e-01, 4.5897200000e-01, 2.4492300000e-01},
      doubles_t{4.9330000000e+01, 1.1370000000e+01, 3.4350000000e+00,
                1.1820000000e+00, 4.1730000000e-01, 1.4280000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1730000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4280000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1167000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6930000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8370000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.6800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.9720000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0270000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.8500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4270000000e+00}));
    return abs_t(name, 7, r0, shells.begin(), shells.end());
} // cc_dash_pwcvqz_7

} // namespace chemcache
