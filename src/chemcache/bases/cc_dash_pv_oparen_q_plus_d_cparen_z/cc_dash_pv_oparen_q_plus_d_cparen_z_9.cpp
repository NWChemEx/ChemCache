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

#include "cc_dash_pv_oparen_q_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_q_plus_d_cparen_z_9() {
    // Basis Set name and origin point
    std::string name("cc-pv(q+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.5000000000e-05, 7.3800000000e-04, 3.8580000000e-03,
                1.5926000000e-02, 5.4289000000e-02, 1.4951300000e-01,
                3.0825200000e-01, 3.9485300000e-01, 2.1103100000e-01,
                1.7151000000e-02, -2.0150000000e-03, 8.6900000000e-04},
      doubles_t{7.4530000000e+04, 1.1170000000e+04, 2.5430000000e+03,
                7.2100000000e+02, 2.3590000000e+02, 8.5600000000e+01,
                3.3550000000e+01, 1.3930000000e+01, 5.9150000000e+00,
                1.8430000000e+00, 7.1240000000e-01, 2.6370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8430000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.1240000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.2000000000e-05, -1.7200000000e-04, -8.9100000000e-04,
                -3.7480000000e-03, -1.2862000000e-02, -3.8061000000e-02,
                -8.6239000000e-02, -1.5586500000e-01, -1.1091400000e-01,
                2.9876100000e-01, 5.8501300000e-01, 2.7115900000e-01},
      doubles_t{7.4530000000e+04, 1.1170000000e+04, 2.5430000000e+03,
                7.2100000000e+02, 2.3590000000e+02, 8.5600000000e+01,
                3.3550000000e+01, 1.3930000000e+01, 5.9150000000e+00,
                1.8430000000e+00, 7.1240000000e-01, 2.6370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6370000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9530000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{6.3470000000e-03, 4.4204000000e-02, 1.6851400000e-01,
                3.6156300000e-01, 4.4217800000e-01, 2.4343500000e-01},
      doubles_t{8.0390000000e+01, 1.8630000000e+01, 5.6940000000e+00,
                1.9530000000e+00, 6.7020000000e-01, 2.1660000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7020000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1660000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0140000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7250000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.5620000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1480000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3760000000e+00}));
    return abs_t(name, 9, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_q_plus_d_cparen_z_9

} // namespace chemcache
