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

#include "cc_dash_pwcvdz_dash_rifit.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pwcvdz_dash_rifit_35() {
    // Basis Set name and origin point
    std::string name("cc-pwcvdz-rifit");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1756500000e+03}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0093800000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3372946380e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2558500000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.8694892430e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8102417460e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2528348080e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.4483300000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0398266900e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0301917400e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{8.9705633000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4806404000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6874940000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.3624200000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9266500000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0337895173e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5714000000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3758197618e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.7568960248e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3387300000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4043846801e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6029231208e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7529190587e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4308210759e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1198329377e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9602700000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3130713064e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{8.3800500000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.6494491209e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0284500000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.4020341079e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.3476400000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0067841591e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4367392054e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.4036289088e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0605877594e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0041400000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1160100000e+02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5860700000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2175100000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4377400000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9902200000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.0491200000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5699500000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6503600000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{6.3972400000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2694000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4457700000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9943400000e+01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{9.8886400000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1018400000e+00}));
    return abs_t(name, 35, r0, shells.begin(), shells.end());
} // cc_dash_pwcvdz_dash_rifit_35

} // namespace chemcache
