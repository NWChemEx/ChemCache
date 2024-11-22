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

#include "cc_dash_pcvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pcvtz_11() {
    // Basis Set name and origin point
    std::string name("cc-pcvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{1.8061800000e-05, 1.4043000000e-04, 7.3843800000e-04,
                           3.1118200000e-03, 1.1208100000e-02, 3.5282800000e-02,
                           9.5989700000e-02, 2.1373500000e-01, 3.4868800000e-01,
                           3.2456600000e-01, 1.1263300000e-01, 7.0679700000e-03,
                           5.9801000000e-04, 2.0413000000e-04,
                           -5.3087000000e-06, 1.8260500000e-06},
                 doubles_t{4.2300000000e+05, 6.3340000000e+04, 1.4410000000e+04,
                           4.0770000000e+03, 1.3280000000e+03, 4.7860000000e+02,
                           1.8620000000e+02, 7.6920000000e+01, 3.3320000000e+01,
                           1.5000000000e+01, 6.8690000000e+00, 2.6830000000e+00,
                           1.1090000000e+00, 4.5400000000e-01, 6.0150000000e-02,
                           2.3820000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.4065300000e-06, -3.4344300000e-05, -1.8011400000e-04,
                -7.6390000000e-04, -2.7524800000e-03, -8.8601600000e-03,
                -2.4793900000e-02, -6.0599500000e-02, -1.1644600000e-01,
                -1.6243700000e-01, -4.3889100000e-02, 3.3791700000e-01,
                5.6134700000e-01, 2.4524000000e-01, 4.0675400000e-03,
                -1.2674600000e-03},
      doubles_t{4.2300000000e+05, 6.3340000000e+04, 1.4410000000e+04,
                4.0770000000e+03, 1.3280000000e+03, 4.7860000000e+02,
                1.8620000000e+02, 7.6920000000e+01, 3.3320000000e+01,
                1.5000000000e+01, 6.8690000000e+00, 2.6830000000e+00,
                1.1090000000e+00, 4.5400000000e-01, 6.0150000000e-02,
                2.3820000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5400000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{6.6301900000e-07, 5.1576900000e-06, 2.7125000000e-05,
                1.1463500000e-04, 4.1511800000e-04, 1.3297800000e-03,
                3.7559500000e-03, 9.1402500000e-03, 1.7985900000e-02,
                2.5147700000e-02, 7.6352200000e-03, -6.1458900000e-02,
                -1.1572100000e-01, -1.5890600000e-01, 6.2640600000e-01,
                4.7540200000e-01},
      doubles_t{4.2300000000e+05, 6.3340000000e+04, 1.4410000000e+04,
                4.0770000000e+03, 1.3280000000e+03, 4.7860000000e+02,
                1.8620000000e+02, 7.6920000000e+01, 3.3320000000e+01,
                1.5000000000e+01, 6.8690000000e+00, 2.6830000000e+00,
                1.1090000000e+00, 4.5400000000e-01, 6.0150000000e-02,
                2.3820000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.3820000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.1890000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2600000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5690000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{2.2439200000e-03, 1.7399700000e-02, 7.7412500000e-02,
                           2.1910200000e-01, 3.7852200000e-01, 3.9490200000e-01,
                           1.6042400000e-01, 2.3331100000e-03, 1.9953600000e-03,
                           -7.7734400000e-04},
                 doubles_t{2.4330000000e+02, 5.7390000000e+01, 1.8100000000e+01,
                           6.5750000000e+00, 2.5210000000e+00, 9.6070000000e-01,
                           3.5120000000e-01, 9.8270000000e-02, 3.7340000000e-02,
                           1.5000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7340000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-2.2240100000e-04, -1.7427700000e-03, -7.7545600000e-03,
                -2.2518700000e-02, -3.8433000000e-02, -4.5017700000e-02,
                -1.9213200000e-02, 1.8269700000e-01, 5.5789700000e-01,
                3.7302200000e-01},
      doubles_t{2.4330000000e+02, 5.7390000000e+01, 1.8100000000e+01,
                6.5750000000e+00, 2.5210000000e+00, 9.6070000000e-01,
                3.5120000000e-01, 9.8270000000e-02, 3.7340000000e-02,
                1.5000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.4040000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5300000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3670000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{6.3600000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4650000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.3970000000e-01}));
    return abs_t(name, 11, r0, shells.begin(), shells.end());
} // cc_dash_pcvtz_11

} // namespace chemcache
