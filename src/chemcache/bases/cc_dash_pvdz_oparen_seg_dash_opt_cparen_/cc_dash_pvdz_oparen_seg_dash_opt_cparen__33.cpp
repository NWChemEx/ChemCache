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

#include "cc_dash_pvdz_oparen_seg_dash_opt_cparen_.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvdz_oparen_seg_dash_opt_cparen__33() {
    // Basis Set name and origin point
    std::string name("cc-pvdz(seg-opt)");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{2.0238000000e-04, 1.5705900000e-03, 8.1662100000e-03,
                           3.3339900000e-02, 1.0926362000e-01, 2.7242387000e-01,
                           4.2391339000e-01, 2.8272958000e-01, 1.1791250000e-02,
                           -3.1083360000e-02},
                 doubles_t{5.5958379000e+05, 8.3879330000e+04, 1.9092668000e+04,
                           5.4073925000e+03, 1.7637559000e+03, 6.3645672000e+02,
                           2.4808843000e+02, 1.0157851000e+02, 3.1475513000e+01,
                           1.3437282000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.8000000000e-07, 8.3100000000e-06, -5.8000000000e-07,
                -2.4182700000e-03, -1.4256540000e-02, -8.1144620000e-02,
                -9.2510290000e-02, 4.7988082000e-01, 6.3084716000e-01,
                8.6917760000e-02},
      doubles_t{5.5958379000e+05, 8.3879330000e+04, 1.9092668000e+04,
                1.7637559000e+03, 6.3645672000e+02, 2.4808843000e+02,
                1.0157851000e+02, 3.1475513000e+01, 1.3437282000e+01,
                4.0086900000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9000000000e-07, -1.1000000000e-06, -9.9400000000e-06,
                1.3484400000e-03, 5.6417400000e-03, 1.4936880000e-02,
                -1.2061760000e-01, -3.1326128000e-01, 5.7263060000e-01,
                7.2048570000e-01},
      doubles_t{5.5958379000e+05, 8.3879330000e+04, 1.9092668000e+04,
                6.3645672000e+02, 2.4808843000e+02, 1.0157851000e+02,
                3.1475513000e+01, 1.3437282000e+01, 4.0086900000e+00,
                1.6849290000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.0000000000e-08, -3.0000000000e-08, -3.0000000000e-08,
                2.8430000000e-05, 8.7250000000e-05, -8.8630000000e-04,
                -1.7324000000e-03, -1.9193300000e-02, -2.0388461000e-01,
                6.7269177000e-01},
      doubles_t{5.5958379000e+05, 8.3879330000e+04, 1.9092668000e+04,
                2.4808843000e+02, 1.0157851000e+02, 3.1475513000e+01,
                1.3437282000e+01, 4.0086900000e+00, 1.6849290000e+00,
                3.0001900000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1358700000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{1.4097000000e-03, 1.1827830000e-02, 5.9332250000e-02,
                           1.9654567000e-01, 3.9993798000e-01, 4.0055899000e-01,
                           1.1485915000e-01, -2.9199500000e-03},
                 doubles_t{3.8863564000e+03, 9.2120201000e+02, 2.9719319000e+02,
                           1.1197508000e+02, 4.6034621000e+01, 1.9874194000e+01,
                           8.3860880000e+00, 3.5587280000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-2.4430000000e-05, -7.1193000000e-04, -5.7848100000e-03,
                -2.4919120000e-02, -1.6075030000e-02, 2.6198036000e-01,
                6.0175036000e-01, 3.6486889000e-01},
      doubles_t{9.2120201000e+02, 2.9719319000e+02, 1.1197508000e+02,
                4.6034621000e+01, 1.9874194000e+01, 8.3860880000e+00,
                3.5587280000e+00, 1.4472820000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.4703000000e-04, 4.0161400000e-03, 1.9430800000e-02,
                9.5619800000e-03, -1.9816200000e-01, -4.8996142000e-01,
                -2.6335171000e-01, 5.2492249000e-01},
      doubles_t{2.9719319000e+02, 1.1197508000e+02, 4.6034621000e+01,
                1.9874194000e+01, 8.3860880000e+00, 3.5587280000e+00,
                1.4472820000e+00, 3.4777900000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0769900000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.4528800000e-02, 1.4113400000e-01, 3.6875790000e-01,
                           4.8406260000e-01, 2.8244340000e-01},
                 doubles_t{8.4424234000e+01, 2.4181589000e+01, 8.4017770000e+00,
                           2.9805020000e+00, 9.7900300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0980000000e-01}));
    return abs_t(name, 33, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_oparen_seg_dash_opt_cparen__33

} // namespace chemcache
