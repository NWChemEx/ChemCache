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

#include "d_dash_aug_dash_cc_dash_pvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t d_dash_aug_dash_cc_dash_pvqz_8() {
    // Basis Set name and origin point
    std::string name("d-aug-cc-pvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.0000000000e-05, 6.9800000000e-04, 3.6640000000e-03,
                1.5218000000e-02, 5.2423000000e-02, 1.4592100000e-01,
                3.0525800000e-01, 3.9850800000e-01, 2.1698000000e-01,
                1.7594000000e-02, -2.5020000000e-03, 9.5400000000e-04},
      doubles_t{6.1420000000e+04, 9.1990000000e+03, 2.0910000000e+03,
                5.9090000000e+02, 1.9230000000e+02, 6.9320000000e+01,
                2.6970000000e+01, 1.1100000000e+01, 4.6820000000e+00,
                1.4280000000e+00, 5.5470000000e-01, 2.0670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4280000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5470000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.0000000000e-05, -1.5900000000e-04, -8.2900000000e-04,
                -3.5080000000e-03, -1.2156000000e-02, -3.6261000000e-02,
                -8.2992000000e-02, -1.5209000000e-01, -1.1533100000e-01,
                2.8897900000e-01, 5.8612800000e-01, 2.7762400000e-01},
      doubles_t{6.1420000000e+04, 9.1990000000e+03, 2.0910000000e+03,
                5.9090000000e+02, 1.9230000000e+02, 6.9320000000e+01,
                2.6970000000e+01, 1.1100000000e+01, 4.6820000000e+00,
                1.4280000000e+00, 5.5470000000e-01, 2.0670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.9590000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5310000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{6.0440000000e-03, 4.1799000000e-02, 1.6114300000e-01,
                3.5673100000e-01, 4.4830900000e-01, 2.4494000000e-01},
      doubles_t{6.3420000000e+01, 1.4660000000e+01, 4.4590000000e+00,
                1.5310000000e+00, 5.3020000000e-01, 1.7500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3020000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3480000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6300000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7750000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3000000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6660000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8460000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1400000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7600000000e-01}));
    return abs_t(name, 8, r0, shells.begin(), shells.end());
} // d_dash_aug_dash_cc_dash_pvqz_8

} // namespace chemcache
