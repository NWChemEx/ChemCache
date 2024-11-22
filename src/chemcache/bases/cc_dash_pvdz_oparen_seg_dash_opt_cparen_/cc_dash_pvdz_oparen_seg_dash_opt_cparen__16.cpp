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

abs_t cc_dash_pvdz_oparen_seg_dash_opt_cparen__16() {
    // Basis Set name and origin point
    std::string name("cc-pvdz(seg-opt)");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4764000000e-04, 1.9203000000e-03, 9.9618800000e-03,
                4.0297500000e-02, 1.2859135000e-01, 3.0341538000e-01,
                4.2107356000e-01, 2.3051703000e-01, 2.0238940000e-02},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                1.0710000000e+03, 3.4980000000e+02, 1.2630000000e+02,
                4.9260000000e+01, 2.0160000000e+01, 5.7200000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.3000000000e-07, 9.1700000000e-06, -1.1440000000e-05,
                -2.9323600000e-03, -1.4719290000e-02, -8.2133190000e-02,
                -5.9246880000e-02, 5.2240927000e-01, 6.2040624000e-01},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                3.4980000000e+02, 1.2630000000e+02, 4.9260000000e+01,
                2.0160000000e+01, 5.7200000000e+00, 2.1820000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.7000000000e-07, -8.0000000000e-07, -6.9800000000e-06,
                8.7289000000e-04, 3.2866800000e-03, 6.4719000000e-03,
                -6.7925580000e-02, -2.2766568000e-01, 6.4363158000e-01},
      doubles_t{1.1080000000e+05, 1.6610000000e+04, 3.7810000000e+03,
                1.2630000000e+02, 4.9260000000e+01, 2.0160000000e+01,
                5.7200000000e+00, 2.1820000000e+00, 4.3270000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5700000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.4754100000e-03, 3.4166830000e-02, 1.4427929000e-01,
                3.5395563000e-01, 4.5960531000e-01, 2.0486597000e-01},
      doubles_t{3.9970000000e+02, 9.4190000000e+01, 2.9750000000e+01,
                1.0770000000e+01, 4.1190000000e+00, 1.6250000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.1940000000e-04, -1.6189500000e-03, -1.5279400000e-03,
                -2.8744400000e-02, 8.3799430000e-02, 5.6422629000e-01},
      doubles_t{9.4190000000e+01, 2.9750000000e+01, 1.0770000000e+01,
                4.1190000000e+00, 1.6250000000e+00, 4.7260000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4070000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.7900000000e-01}));
    return abs_t(name, 16, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_oparen_seg_dash_opt_cparen__16

} // namespace chemcache
