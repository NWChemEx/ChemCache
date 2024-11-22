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

abs_t mini_19() {
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
      doubles_t{6.4874696171e-02, 3.8085927752e-01, 6.7736806002e-01},
      doubles_t{1.7211755000e+03, 2.6001633000e+02, 5.6624554000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-9.7873101797e-02, 6.5955981211e-01, 4.0652950746e-01},
      doubles_t{7.5055600000e+01, 6.6911626000e+00, 2.6671665000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.1293020289e-01, 6.8924910934e-01, 4.3392540588e-01},
      doubles_t{4.6689367000e+00, 7.0001330000e-01, 2.7533360000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.5298000020e-01, 6.8172840089e-01, 4.1448290054e-01},
      doubles_t{2.5231150000e-01, 3.7631100000e-02, 1.6218300000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{9.8035402529e-02, 4.5104081163e-01, 6.1770761593e-01},
      doubles_t{8.5789846000e+01, 1.9254794000e+01, 5.2686239000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.3863598893e-01, 5.7105677351e-01, 3.1653568531e-01},
      doubles_t{1.6831452000e+00, 6.2580940000e-01, 2.2389800000e-01}));
    return abs_t(name, 19, r0, shells.begin(), shells.end());
} // mini_19

} // namespace chemcache
