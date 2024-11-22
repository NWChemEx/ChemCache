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

abs_t mini_67() {
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
      doubles_t{6.0408602586e-02, 3.6425931559e-01, 6.9341642968e-01},
      doubles_t{2.3180533000e+04, 3.5116702000e+03, 7.7332413000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1168588873e-01, 7.3114462620e-01, 3.3544556614e-01},
      doubles_t{1.0292619000e+03, 1.0170128000e+02, 4.3328224000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9150659807e-01, 8.7030109425e-01, 2.9473209805e-01},
      doubles_t{9.2657169000e+01, 1.7937269000e+01, 8.4777760000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.7799388358e-01, 7.1887096876e-01, 5.0973727785e-01},
      doubles_t{1.6692436000e+01, 4.3358910000e+00, 2.2322713000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.7470994320e-03, 2.2930948069e-01, 7.8843143362e-01},
      doubles_t{4.8851960000e+00, 7.9382900000e-01, 3.9985840000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{7.9827304288e-02, 2.5017381344e-01, 7.5698434067e-01},
      doubles_t{9.4909290000e-01, 7.3420400000e-02, 2.7928600000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{7.7398498691e-02, 4.1554409297e-01, 6.4676118906e-01},
      doubles_t{1.5737315000e+03, 3.6996920000e+02, 1.1106992000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.5412699260e-02, 3.8895769187e-01, 6.7788848583e-01},
      doubles_t{2.0941284000e+02, 4.1838408000e+01, 1.6853075000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.3292881804e-01, 5.4753422966e-01, 1.9345931048e-01},
      doubles_t{8.1634159000e+00, 3.8567519000e+00, 1.8073015000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.6625099344e-02, 5.3835387875e-01, 5.2502677928e-01},
      doubles_t{2.1201059000e+00, 8.2864220000e-01, 2.9950200000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.0852549811e-01, 4.7707619170e-01, 5.9045508973e-01},
      doubles_t{2.3721490000e+02, 6.7889543000e+01, 2.2617285000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.8585568511e-01, 5.8195536969e-01, 2.9478708465e-01},
      doubles_t{1.3093652000e+01, 4.9213964000e+00, 1.8661708000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 3,
                 doubles_t{8.3792200063e-02, 3.2831860025e-01, 5.2332550039e-01,
                           4.2423540032e-01},
                 doubles_t{4.4012473000e+01, 1.3416270000e+01, 4.4270768000e+00,
                           1.2896098000e+00}));
    return abs_t(name, 67, r0, shells.begin(), shells.end());
} // mini_67

} // namespace chemcache
