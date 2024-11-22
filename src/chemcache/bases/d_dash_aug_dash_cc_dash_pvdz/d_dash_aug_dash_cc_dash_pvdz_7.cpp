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

#include "d_dash_aug_dash_cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t d_dash_aug_dash_cc_dash_pvdz_7() {
    // Basis Set name and origin point
    std::string name("d-aug-cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.0000000000e-04, 5.3890000000e-03, 2.7406000000e-02,
                1.0320700000e-01, 2.7872300000e-01, 4.4854000000e-01,
                2.7823800000e-01, 1.5440000000e-02, -2.8640000000e-03},
      doubles_t{9.0460000000e+03, 1.3570000000e+03, 3.0930000000e+02,
                8.7730000000e+01, 2.8560000000e+01, 1.0210000000e+01,
                3.8380000000e+00, 7.4660000000e-01, 2.2480000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.5300000000e-04, -1.2080000000e-03, -5.9920000000e-03,
                -2.4544000000e-02, -6.7459000000e-02, -1.5807800000e-01,
                -1.2183100000e-01, 5.4900300000e-01, 5.7881500000e-01},
      doubles_t{9.0460000000e+03, 1.3570000000e+03, 3.0930000000e+02,
                8.7730000000e+01, 2.8560000000e+01, 1.0210000000e+01,
                3.8380000000e+00, 7.4660000000e-01, 2.2480000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2480000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1240000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6700000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{3.9919000000e-02, 2.1716900000e-01, 5.1031900000e-01,
                           4.6221400000e-01},
                 doubles_t{1.3550000000e+01, 2.9170000000e+00, 7.9730000000e-01,
                           2.1850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1850000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.6110000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.1700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4700000000e-02}));
    return abs_t(name, 7, r0, shells.begin(), shells.end());
} // d_dash_aug_dash_cc_dash_pvdz_7

} // namespace chemcache
