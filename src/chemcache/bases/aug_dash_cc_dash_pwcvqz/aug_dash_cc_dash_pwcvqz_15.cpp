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

#include "aug_dash_cc_dash_pwcvqz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pwcvqz_15() {
    // Basis Set name and origin point
    std::string name("aug-cc-pwcvqz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4745000000e-05, 1.9246500000e-04, 1.0120200000e-03,
                4.2726100000e-03, 1.5416100000e-02, 4.8597600000e-02,
                1.3006000000e-01, 2.7451400000e-01, 3.8540200000e-01,
                2.5593400000e-01, 3.9123700000e-02, -3.6801000000e-03,
                2.0821100000e-03, -7.8847900000e-04, 4.5405400000e-04,
                -1.2664000000e-04},
      doubles_t{6.1520000000e+05, 9.2120000000e+04, 2.0950000000e+04,
                5.9200000000e+03, 1.9220000000e+03, 6.8800000000e+02,
                2.6500000000e+02, 1.0820000000e+02, 4.6220000000e+01,
                2.0230000000e+01, 7.8590000000e+00, 3.5470000000e+00,
                1.5640000000e+00, 4.8880000000e-01, 2.2660000000e-01,
                9.3310000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.7220500000e-06, -5.2231100000e-05, -2.7536100000e-04,
                -1.1630700000e-03, -4.2428100000e-03, -1.3611400000e-02,
                -3.8511400000e-02, -9.0664300000e-02, -1.6658400000e-01,
                -1.6144700000e-01, 1.4678100000e-01, 5.6668200000e-01,
                4.1643300000e-01, 3.4384400000e-02, -7.8063800000e-03,
                1.9225900000e-03},
      doubles_t{6.1520000000e+05, 9.2120000000e+04, 2.0950000000e+04,
                5.9200000000e+03, 1.9220000000e+03, 6.8800000000e+02,
                2.6500000000e+02, 1.0820000000e+02, 4.6220000000e+01,
                2.0230000000e+01, 7.8590000000e+00, 3.5470000000e+00,
                1.5640000000e+00, 4.8880000000e-01, 2.2660000000e-01,
                9.3310000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8880000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2660000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.8474000000e-06, 1.4338000000e-05, 7.5722800000e-05,
                3.1920500000e-04, 1.1685100000e-03, 3.7426700000e-03,
                1.0681700000e-02, 2.5265700000e-02, 4.7928300000e-02,
                4.7709600000e-02, -4.6652500000e-02, -2.3496800000e-01,
                -3.1133700000e-01, 2.5710900000e-01, 6.5365500000e-01,
                2.9421200000e-01},
      doubles_t{6.1520000000e+05, 9.2120000000e+04, 2.0950000000e+04,
                5.9200000000e+03, 1.9220000000e+03, 6.8800000000e+02,
                2.6500000000e+02, 1.0820000000e+02, 4.6220000000e+01,
                2.0230000000e+01, 7.8590000000e+00, 3.5470000000e+00,
                1.5640000000e+00, 4.8880000000e-01, 2.2660000000e-01,
                9.3310000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.3310000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3298000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4340000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1130000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4060000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.6150000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{4.2101500000e-04, 3.6098500000e-03, 1.8921700000e-02,
                           7.0556000000e-02, 1.8815700000e-01, 3.3870900000e-01,
                           3.8194300000e-01, 1.9526100000e-01, 1.9941600000e-02,
                           -1.3512100000e-03, 5.1714400000e-04},
                 doubles_t{1.3670000000e+03, 3.2400000000e+02, 1.0460000000e+02,
                           3.9370000000e+01, 1.6260000000e+01, 7.0560000000e+00,
                           3.1300000000e+00, 1.3940000000e+00, 5.1790000000e-01,
                           2.0320000000e-01, 7.6980000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1790000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0320000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.0082700000e-04, -8.5449900000e-04, -4.5711600000e-03,
                -1.7032700000e-02, -4.7520400000e-02, -8.5278600000e-02,
                -1.0967600000e-01, -1.6118100000e-02, 3.2289300000e-01,
                5.4573800000e-01, 2.6853800000e-01},
      doubles_t{1.3670000000e+03, 3.2400000000e+02, 1.0460000000e+02,
                3.9370000000e+01, 1.6260000000e+01, 7.0560000000e+00,
                3.1300000000e+00, 1.3940000000e+00, 5.1790000000e-01,
                2.0320000000e-01, 7.6980000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.6980000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1120000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7200000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9566000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1190000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5900000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0360000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9400000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0518000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3100000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0300000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8000000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0900000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3280000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5000000000e-01}));
    return abs_t(name, 15, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pwcvqz_15

} // namespace chemcache
