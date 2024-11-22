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

#include "cc_dash_pcv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pcv5z_20() {
    // Basis Set name and origin point
    std::string name("cc-pcv5z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{4.3000000000e-07, 3.2900000000e-06, 1.7140000000e-05,
                           7.2280000000e-05, 2.6522000000e-04, 8.7484000000e-04,
                           2.6499100000e-03, 7.4844700000e-03, 1.9730930000e-02,
                           4.7797990000e-02, 1.0359562000e-01, 1.9293018000e-01,
                           2.8626757000e-01, 2.9278494000e-01, 1.6312798000e-01,
                           3.2216880000e-02, 5.9258000000e-04},
                 doubles_t{2.8249600000e+07, 4.2501900000e+06, 9.7501400000e+05,
                           2.7744600000e+05, 9.0454700000e+04, 3.2515900000e+04,
                           1.2610600000e+04, 5.1907700000e+03, 2.2423300000e+03,
                           1.0115300000e+03, 4.7554700000e+02, 2.3207200000e+02,
                           1.1687400000e+02, 6.0342700000e+01, 3.1524400000e+01,
                           1.6030800000e+01, 8.4599200000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.2000000000e-07, -9.6000000000e-07, -4.9700000000e-06,
                -2.0990000000e-05, -7.7010000000e-05, -2.5431000000e-04,
                -7.7168000000e-04, -2.1908300000e-03, -5.8356700000e-03,
                -1.4468690000e-02, -3.2803840000e-02, -6.6744350000e-02,
                -1.1668667000e-01, -1.5919701000e-01, -1.1620650000e-01,
                1.0547012000e-01, 4.2483472000e-01},
      doubles_t{2.8249600000e+07, 4.2501900000e+06, 9.7501400000e+05,
                2.7744600000e+05, 9.0454700000e+04, 3.2515900000e+04,
                1.2610600000e+04, 5.1907700000e+03, 2.2423300000e+03,
                1.0115300000e+03, 4.7554700000e+02, 2.3207200000e+02,
                1.1687400000e+02, 6.0342700000e+01, 3.1524400000e+01,
                1.6030800000e+01, 8.4599200000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 0,
                 doubles_t{4.0000000000e-08, 3.3000000000e-07, 1.7200000000e-06,
                           7.2400000000e-06, 2.6560000000e-05, 8.7720000000e-05,
                           2.6618000000e-04, 7.5622000000e-04, 2.0161800000e-03,
                           5.0120300000e-03, 1.1414850000e-02, 2.3454260000e-02,
                           4.1737840000e-02, 5.8943050000e-02, 4.5034190000e-02,
                           -4.4376180000e-02, -2.1875449000e-01},
                 doubles_t{2.8249600000e+07, 4.2501900000e+06, 9.7501400000e+05,
                           2.7744600000e+05, 9.0454700000e+04, 3.2515900000e+04,
                           1.2610600000e+04, 5.1907700000e+03, 2.2423300000e+03,
                           1.0115300000e+03, 4.7554700000e+02, 2.3207200000e+02,
                           1.1687400000e+02, 6.0342700000e+01, 3.1524400000e+01,
                           1.6030800000e+01, 8.4599200000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5271800000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4185900000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2512900000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{6.4343000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.2918000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4239000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{7.8190000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.7630000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.8290000000e-02}));
    shells.emplace_back(
      make_shell(pure_t::pure, 1,
                 doubles_t{2.3000000000e-07, 2.3930000000e-05, 2.1218000000e-04,
                           1.2273300000e-03, 5.4527700000e-03, 1.9626690000e-02,
                           5.8279930000e-02, 1.3955299000e-01, 2.5592562000e-01,
                           3.3998132000e-01},
                 doubles_t{4.0635300000e+05, 1.3600700000e+04, 3.2354700000e+03,
                           1.0530500000e+03, 4.0356900000e+02, 1.7149000000e+02,
                           7.8236300000e+01, 3.7606500000e+01, 1.8709000000e+01,
                           9.5004100000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-7.0000000000e-08, -7.7900000000e-06, -6.9080000000e-05,
                -4.0059000000e-04, -1.7849000000e-03, -6.4832000000e-03,
                -1.9562440000e-02, -4.8247670000e-02, -9.1286110000e-02,
                -1.2842086000e-01},
      doubles_t{4.0635300000e+05, 1.3600700000e+04, 3.2354700000e+03,
                1.0530500000e+03, 4.0356900000e+02, 1.7149000000e+02,
                7.8236300000e+01, 3.7606500000e+01, 1.8709000000e+01,
                9.5004100000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.9055700000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5267000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2713900000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.2601000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0099000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1769000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1160000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1420000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.9047000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1914000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.4460000000e-03, 2.1360000000e-02, 7.5448000000e-02,
                           1.7152800000e-01, 2.6061900000e-01, 3.3474500000e-01,
                           3.7814300000e-01, 2.1650400000e-01},
                 doubles_t{3.9562800000e+01, 1.1437300000e+01, 3.9674300000e+00,
                           1.5247800000e+00, 5.9047000000e-01, 2.1914000000e-01,
                           7.9570000000e-02, 2.8340000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9570000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.8340000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.7410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5800000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6000000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{3.0230000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0120000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5340000000e-01}));
    return abs_t(name, 20, r0, shells.begin(), shells.end());
} // cc_dash_pcv5z_20

} // namespace chemcache
