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

#include "cc_dash_pvtz_dash_rifit.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvtz_dash_rifit_14() {
    // Basis Set name and origin point
    std::string name("cc-pvtz-rifit");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5489886368e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{9.3369336022e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4848188316e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1101816659e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9645817082e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7827045121e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1314513609e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2521096076e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8105221542e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6363612894e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.3831114156e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8657375575e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7885088086e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9711361601e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1225083205e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5040198463e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3335561695e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5281212451e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3801820114e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.6761446485e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6107939018e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9368390161e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1948646470e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.6118405469e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0514700903e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9401974957e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5489998282e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3000060620e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.1563942171e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9155067749e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8443341045e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1224369309e-01}));
    return abs_t(name, 14, r0, shells.begin(), shells.end());
} // cc_dash_pvtz_dash_rifit_14

} // namespace chemcache
