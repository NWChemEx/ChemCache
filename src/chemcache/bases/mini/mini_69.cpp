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

#include "mini.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t mini_69() {
    // Basis Set name and origin point
    std::string name("mini");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.0240802580e-02, 3.6383741558e-01, 6.9391762973e-01},
      doubles_t{2.4656410000e+04, 3.7325062000e+03, 8.2163815000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1203019610e-01, 7.2673947469e-01, 3.3976558817e-01},
      doubles_t{1.0914725000e+03, 1.0854804000e+02, 4.6629354000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9244747394e-01, 8.7251922226e-01, 2.9461427375e-01},
      doubles_t{9.8840728000e+01, 1.9274990000e+01, 8.9158303000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.7863799745e-01, 7.0195659528e-01, 5.2820799645e-01},
      doubles_t{1.7990337000e+01, 4.7415251000e+00, 2.4155100000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.2713000977e-02, 2.4587301889e-01, 7.7580765960e-01},
      doubles_t{4.9925787000e+00, 8.3566170000e-01, 4.1594390000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.9699499802e-02, 2.2358509945e-01, 7.8195759806e-01},
      doubles_t{1.0216999000e+00, 7.7869100000e-02, 2.9620900000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{7.6625504082e-02, 4.1378282204e-01, 6.4874443456e-01},
      doubles_t{1.6854211000e+03, 3.9599853000e+02, 1.1891754000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.5908497815e-02, 3.9232147612e-01, 6.7402855898e-01},
      doubles_t{2.2782181000e+02, 4.4528319000e+01, 1.8046201000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.2021347871e-01, 5.5337056321e-01, 2.0161598660e-01},
      doubles_t{8.9691876000e+00, 4.2464848000e+00, 1.9503073000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.1915899575e-02, 5.4211328948e-01, 5.1837198994e-01},
      doubles_t{2.1349230000e+00, 8.6763440000e-01, 3.1114730000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.0686660036e-01, 4.7574480161e-01, 5.9238660201e-01},
      doubles_t{2.5688245000e+02, 7.3519651000e+01, 2.4520879000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.8935690966e-01, 5.7839341931e-01, 2.9463840984e-01},
      doubles_t{1.4183386000e+01, 5.3348920000e+00, 2.0382571000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3,
                 doubles_t{8.8021197376e-02, 3.3680468996e-01, 5.2277058442e-01,
                           4.1597618760e-01},
                 doubles_t{4.6858648000e+01, 1.4376483000e+01, 4.7433937000e+00,
                           1.3719529000e+00}));
    return abs_t(name, 69, r0, shells.begin(), shells.end());
} // mini_69

} // namespace chemcache
