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

abs_t mini_70() {
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
      doubles_t{6.0217797372e-02, 3.6372048413e-01, 6.9402446971e-01},
      doubles_t{2.5385694000e+04, 3.8436032000e+03, 8.4615684000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1201189835e-01, 7.3083348922e-01, 3.3650889503e-01},
      doubles_t{1.1234055000e+03, 1.1169297000e+02, 4.7233256000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9465308853e-01, 8.8516946554e-01, 2.8258658900e-01},
      doubles_t{1.0178796000e+02, 1.9813955000e+01, 9.1947339000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.7032011500e-01, 7.2190452924e-01, 5.0080532029e-01},
      doubles_t{1.8855724000e+01, 4.8006027000e+00, 2.4675257000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.9842199412e-02, 2.1967679349e-01, 8.0312317619e-01},
      doubles_t{5.1383398000e+00, 8.6005930000e-01, 4.3522800000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.4708803141e-02, 2.0341470855e-01, 8.0212113373e-01},
      doubles_t{1.0515290000e+00, 8.0142600000e-02, 3.0485600000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{7.6544501139e-02, 4.1356810615e-01, 6.4890450965e-01},
      doubles_t{1.7377854000e+03, 4.0850276000e+02, 1.2273231000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.6068402063e-02, 3.9315502248e-01, 6.7324723850e-01},
      doubles_t{2.3449301000e+02, 4.6044797000e+01, 1.8679254000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.3541433064e-01, 5.4180414949e-01, 1.9765891806e-01},
      doubles_t{9.1334438000e+00, 4.3362256000e+00, 2.0173973000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.6651400355e-02, 5.4655930728e-01, 5.1110400681e-01},
      doubles_t{2.1972533000e+00, 8.7979160000e-01, 3.1464740000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.0699449588e-01, 4.7611868167e-01, 5.9150457723e-01},
      doubles_t{2.6549751000e+02, 7.6107209000e+01, 2.5437233000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.9220651440e-01, 5.8117412864e-01, 2.8978801428e-01},
      doubles_t{1.4703007000e+01, 5.5153232000e+00, 2.0915122000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3,
                 doubles_t{9.0153299961e-02, 3.4030879985e-01, 5.2107449978e-01,
                           4.1283349982e-01},
                 doubles_t{4.8314683000e+01, 1.4868502000e+01, 4.9219249000e+00,
                           1.4244136000e+00}));
    return abs_t(name, 70, r0, shells.begin(), shells.end());
} // mini_70

} // namespace chemcache
