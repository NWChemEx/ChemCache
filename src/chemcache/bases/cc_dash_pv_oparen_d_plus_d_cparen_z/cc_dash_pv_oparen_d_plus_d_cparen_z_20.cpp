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

#include "cc_dash_pv_oparen_d_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pv_oparen_d_plus_d_cparen_z_20() {
    // Basis Set name and origin point
    std::string name("cc-pv(d+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.2145000000e-04, 1.7183000000e-03, 8.9234800000e-03,
                3.6301830000e-02, 1.1762220000e-01, 2.8604350000e-01,
                4.2260710000e-01, 2.5774370000e-01, 2.3918930000e-02,
                -4.9521800000e-03, 1.7177900000e-03, -8.9209000000e-04,
                2.4510000000e-04, -1.2395000000e-04},
      doubles_t{1.9000070000e+05, 2.8481460000e+04, 6.4827010000e+03,
                1.8358910000e+03, 5.9872430000e+02, 2.1588410000e+02,
                8.4012420000e+01, 3.4224880000e+01, 1.0024970000e+01,
                4.0559200000e+00, 1.0202610000e+00, 4.2686500000e-01,
                6.3347000000e-02, 2.6301000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.4530000000e-05, -4.9662000000e-04, -2.6282600000e-03,
                -1.0668450000e-02, -3.7135090000e-02, -9.8042840000e-02,
                -2.0342690000e-01, -1.5244650000e-01, 4.8279410000e-01,
                6.2923840000e-01, 6.1648420000e-02, -1.4799710000e-02,
                3.6108900000e-03, -1.7927300000e-03},
      doubles_t{1.9000070000e+05, 2.8481460000e+04, 6.4827010000e+03,
                1.8358910000e+03, 5.9872430000e+02, 2.1588410000e+02,
                8.4012420000e+01, 3.4224880000e+01, 1.0024970000e+01,
                4.0559200000e+00, 1.0202610000e+00, 4.2686500000e-01,
                6.3347000000e-02, 2.6301000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.2230000000e-05, 1.7170000000e-04, 9.0452000000e-04,
                3.7034300000e-03, 1.2837500000e-02, 3.4754590000e-02,
                7.3034910000e-02, 6.1000830000e-02, -2.4292930000e-01,
                -4.8708500000e-01, 5.6502800000e-01, 6.5574390000e-01,
                2.6728940000e-02, -9.9995900000e-03},
      doubles_t{1.9000070000e+05, 2.8481460000e+04, 6.4827010000e+03,
                1.8358910000e+03, 5.9872430000e+02, 2.1588410000e+02,
                8.4012420000e+01, 3.4224880000e+01, 1.0024970000e+01,
                4.0559200000e+00, 1.0202610000e+00, 4.2686500000e-01,
                6.3347000000e-02, 2.6301000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.3100000000e-06, 4.1110000000e-05, 2.1568000000e-04,
                8.8827000000e-04, 3.0581300000e-03, 8.3760800000e-03,
                1.7410560000e-02, 1.5154530000e-02, -6.2079190000e-02,
                -1.2611800000e-01, 1.7360690000e-01, 3.7822940000e-01,
                -6.5964700000e-01, -4.9022160000e-01},
      doubles_t{1.9000070000e+05, 2.8481460000e+04, 6.4827010000e+03,
                1.8358910000e+03, 5.9872430000e+02, 2.1588410000e+02,
                8.4012420000e+01, 3.4224880000e+01, 1.0024970000e+01,
                4.0559200000e+00, 1.0202610000e+00, 4.2686500000e-01,
                6.3347000000e-02, 2.6301000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.6301000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.9816600000e-03, 1.6129440000e-02, 7.6578510000e-02,
                2.3269590000e-01, 4.2445210000e-01, 3.7326400000e-01,
                7.8685300000e-02, -5.9992700000e-03, 2.6425700000e-03,
                -8.5694000000e-04, 3.3147000000e-04},
      doubles_t{1.0720430000e+03, 2.5384390000e+02, 8.1316260000e+01,
                3.0241830000e+01, 1.2101100000e+01, 5.0225540000e+00,
                1.9092200000e+00, 7.7130400000e-01, 3.0057000000e-01,
                7.6649000000e-02, 2.7772000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-6.4891000000e-04, -5.2790700000e-03, -2.5811310000e-02,
                -8.0628920000e-02, -1.5846550000e-01, -1.2816820000e-01,
                2.5610100000e-01, 5.8724070000e-01, 3.0372560000e-01,
                1.4164510000e-02, -1.1522400000e-03},
      doubles_t{1.0720430000e+03, 2.5384390000e+02, 8.1316260000e+01,
                3.0241830000e+01, 1.2101100000e+01, 5.0225540000e+00,
                1.9092200000e+00, 7.7130400000e-01, 3.0057000000e-01,
                7.6649000000e-02, 2.7772000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.3595000000e-04, 1.0942000000e-03, 5.4268000000e-03,
                1.6747180000e-02, 3.3898630000e-02, 2.5311830000e-02,
                -5.8957130000e-02, -1.5876120000e-01, -8.5545230000e-02,
                5.4464670000e-01, 5.6631280000e-01},
      doubles_t{1.0720430000e+03, 2.5384390000e+02, 8.1316260000e+01,
                3.0241830000e+01, 1.2101100000e+01, 5.0225540000e+00,
                1.9092200000e+00, 7.7130400000e-01, 3.0057000000e-01,
                7.6649000000e-02, 2.7772000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7772000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.2849000000e-02, 1.4819200000e-01, 3.1092100000e-01,
                           4.5219500000e-01, 4.8086500000e-01},
                 doubles_t{1.0318200000e+01, 2.5924200000e+00, 7.6170000000e-01,
                           2.0838000000e-01, 5.3700000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.3700000000e-02}));
    return abs_t(name, 20, r0, shells.begin(), shells.end());
} // cc_dash_pv_oparen_d_plus_d_cparen_z_20

} // namespace chemcache
