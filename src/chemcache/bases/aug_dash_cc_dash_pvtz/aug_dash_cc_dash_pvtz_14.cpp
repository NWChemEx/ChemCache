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

#include "aug_dash_cc_dash_pvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pvtz_14() {
    // Basis Set name and origin point
    std::string name("aug-cc-pvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.2510100000e-05, 4.8555300000e-04, 2.5451600000e-03,
                1.0586600000e-02, 3.6878700000e-02, 1.0747900000e-01,
                2.4793600000e-01, 3.9092700000e-01, 3.0202600000e-01,
                5.5923600000e-02, -4.0240600000e-03, 2.5803000000e-03,
                -1.0382100000e-03, 6.0793000000e-04, -1.5402200000e-04},
      doubles_t{2.5490000000e+05, 3.8190000000e+04, 8.6900000000e+03,
                2.4620000000e+03, 8.0480000000e+02, 2.9130000000e+02,
                1.1360000000e+02, 4.6750000000e+01, 1.9820000000e+01,
                7.7080000000e+00, 3.3400000000e+00, 1.4020000000e+00,
                4.3870000000e-01, 2.0700000000e-01, 7.9440000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.6637000000e-05, -1.2931000000e-04, -6.7882800000e-04,
                -2.8411700000e-03, -1.0055100000e-02, -3.0577400000e-02,
                -7.7725600000e-02, -1.5423600000e-01, -1.8036800000e-01,
                7.9821800000e-02, 5.4744100000e-01, 4.8011900000e-01,
                4.7484600000e-02, -1.0699600000e-02, 2.1987100000e-03},
      doubles_t{2.5490000000e+05, 3.8190000000e+04, 8.6900000000e+03,
                2.4620000000e+03, 8.0480000000e+02, 2.9130000000e+02,
                1.1360000000e+02, 4.6750000000e+01, 1.9820000000e+01,
                7.7080000000e+00, 3.3400000000e+00, 1.4020000000e+00,
                4.3870000000e-01, 2.0700000000e-01, 7.9440000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3870000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.2625700000e-06, 3.3106200000e-05, 1.7401500000e-04,
                7.2757400000e-04, 2.5833300000e-03, 7.8635400000e-03,
                2.0215500000e-02, 4.0732000000e-02, 4.9935800000e-02,
                -2.4939600000e-02, -1.9035000000e-01, -3.1835000000e-01,
                9.4803600000e-02, 6.8118000000e-01, 3.9567200000e-01},
      doubles_t{2.5490000000e+05, 3.8190000000e+04, 8.6900000000e+03,
                2.4620000000e+03, 8.0480000000e+02, 2.9130000000e+02,
                1.1360000000e+02, 4.6750000000e+01, 1.9820000000e+01,
                7.7080000000e+00, 3.3400000000e+00, 1.4020000000e+00,
                4.3870000000e-01, 2.0700000000e-01, 7.9440000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9440000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3000000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.9204500000e-03, 1.5355200000e-02, 7.1399100000e-02,
                2.1305200000e-01, 3.9035400000e-01, 3.9372100000e-01,
                1.3256500000e-01, 3.9563000000e-03, 3.3162400000e-04},
      doubles_t{4.8150000000e+02, 1.1390000000e+02, 3.6230000000e+01,
                1.3340000000e+01, 5.2520000000e+00, 2.1200000000e+00,
                8.5610000000e-01, 2.5280000000e-01, 7.8890000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.5610000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-4.0522000000e-04, -3.3589600000e-03, -1.5286000000e-02,
                -4.8921800000e-02, -8.5500800000e-02, -1.1213700000e-01,
                6.1827400000e-02, 5.5191900000e-01, 5.2349200000e-01},
      doubles_t{4.8150000000e+02, 1.1390000000e+02, 3.6230000000e+01,
                1.3340000000e+01, 5.2520000000e+00, 2.1200000000e+00,
                8.5610000000e-01, 2.5280000000e-01, 7.8890000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8890000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2500000000e-01}));
    return abs_t(name, 14, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pvtz_14

} // namespace chemcache
