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

#include "seg_dash_cc_dash_pv5z_dash_pp.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t seg_dash_cc_dash_pv5z_dash_pp_77() {
    // Basis Set name and origin point
    std::string name("seg-cc-pv5z-pp");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{4.5955457600e-05, 4.8807639100e-03, -4.6733179200e-02,
                3.5588035500e-01, -1.7284667900e+00, 4.5009590000e+00,
                -2.5571058900e+00},
      doubles_t{3.6116800000e+02, 9.5867400000e+01, 5.9852700000e+01,
                3.6825200000e+01, 2.3013300000e+01, 1.3681000000e+01,
                7.2456100000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 0, doubles_t{9.6181075800e-01, 4.0658874500e-02},
                 doubles_t{4.5281600000e+00, 2.4987800000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 0, doubles_t{4.9604522100e-01, 5.1546588600e-01},
                 doubles_t{1.5530900000e+00, 9.4230200000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4608800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7079300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3475300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4553000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0596000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.0456632600e-04, -1.0141652900e-02, 9.8305689300e-02,
                -4.2904455200e-01, 4.1405595900e-01, 8.8737499300e-01},
      doubles_t{5.3148900000e+01, 3.3180700000e+01, 2.0720300000e+01,
                1.2941100000e+01, 8.0825100000e+00, 4.9702700000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1, doubles_t{4.4237811400e-01, 5.8269415500e-01},
                 doubles_t{1.7771000000e+00, 9.9501000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3654600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7106700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2725200000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8607000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6454000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.9679846700e-04, -2.8108560700e-03, 2.0379330100e-02,
                -1.5376654700e-01, 4.0374568000e-01, 6.8745758000e-01},
      doubles_t{4.3813000000e+01, 1.9899100000e+01, 1.2426200000e+01,
                5.6402700000e+00, 1.8794500000e+00, 1.0692500000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9955400000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2590800000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7070200000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6184000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1037000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8012000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7810000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2810000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0870000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1503000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{9.1390000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8840000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6762000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{6.5620000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 6, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3796000000e+00}));
    return abs_t(name, 77, r0, shells.begin(), shells.end());
} // seg_dash_cc_dash_pv5z_dash_pp_77

} // namespace chemcache
