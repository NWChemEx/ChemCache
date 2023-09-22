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

#include "bases.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
pv6z atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(pv6z_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(pv6z_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("pv6z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(1): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7949240000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.1071600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0480200000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{4.4000000000e-05, 3.7200000000e-04, 2.0940000000e-03,
                        8.8630000000e-03, 3.0540000000e-02, 9.0342000000e-02,
                        2.1323900000e-01, 3.5235000000e-01, 3.3965700000e-01,
                        1.0733000000e-01},
              doubles_t{1.7767755600e+03, 2.5401771200e+02, 5.4698039000e+01,
                        1.5018344000e+01, 4.9150780000e+00, 1.7949240000e+00,
                        7.1071600000e-01, 3.0480200000e-01, 1.3804600000e-01,
                        6.2157000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3804600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.2157000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3950000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9950000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.0600000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1100000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8660000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9740000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2150000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.8600000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3900000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6650000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6190000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.1500000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4480000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3120000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7580000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(6): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{6.0000000000e-06, 4.4000000000e-05, 2.3200000000e-04,
                        9.7900000000e-04, 3.5500000000e-03, 1.1439000000e-02,
                        3.3002000000e-02, 8.3996000000e-02, 1.8068400000e-01,
                        3.0498900000e-01, 3.4128800000e-01},
              doubles_t{3.1222900000e+05, 4.6749690000e+04, 1.0638550000e+04,
                        3.0132870000e+03, 9.8304640000e+02, 3.5488770000e+02,
                        1.3840340000e+02, 5.7361840000e+01, 2.4928380000e+01,
                        1.1229640000e+01, 5.2015490000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.0000000000e-06, 9.0000000000e-06, 4.9000000000e-05,
                        2.0900000000e-04, 7.6000000000e-04, 2.4540000000e-03,
                        7.2020000000e-03, 1.8795000000e-02, 4.3249000000e-02,
                        8.2621000000e-02, 1.2851800000e-01},
              doubles_t{3.1222900000e+05, 4.6749690000e+04, 1.0638550000e+04,
                        3.0132870000e+03, 9.8304640000e+02, 3.5488770000e+02,
                        1.3840340000e+02, 5.7361840000e+01, 2.4928380000e+01,
                        1.1229640000e+01, 5.2015490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4265640000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.6734400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4559900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9712200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.6347000000e-02}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.2000000000e-05, 5.3700000000e-04, 2.9260000000e-03,
                        1.1503000000e-02, 3.5721000000e-02, 9.3120000000e-02},
              doubles_t{4.7864150000e+02, 1.1341840000e+02, 3.6816230000e+01,
                        1.3964920000e+01, 5.7807870000e+00, 2.5435940000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1737790000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5367500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6217300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2361700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.6671000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5584000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3160000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.6500000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0208300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6753500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8629090000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1928790000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9703300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0709700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1260000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3400000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5833300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9055100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9396200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5090000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(7): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{6.0000000000e-06, 4.3000000000e-05, 2.2900000000e-04,
                        9.6500000000e-04, 3.5000000000e-03, 1.1284000000e-02,
                        3.2599000000e-02, 8.3230000000e-02, 1.7993000000e-01,
                        3.0495500000e-01, 3.4115900000e-01},
              doubles_t{4.3238110000e+05, 6.4743420000e+04, 1.4733700000e+04,
                        4.1732540000e+03, 1.3614770000e+03, 4.9150490000e+02,
                        1.9168570000e+02, 7.9454990000e+01, 3.4549680000e+01,
                        1.5586730000e+01, 7.2354440000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.0000000000e-06, 1.0000000000e-05, 5.0000000000e-05,
                        2.1300000000e-04, 7.7500000000e-04, 2.5040000000e-03,
                        7.3620000000e-03, 1.9286000000e-02, 4.4699000000e-02,
                        8.6048000000e-02, 1.3324900000e-01},
              doubles_t{4.3238110000e+05, 6.4743420000e+04, 1.4733700000e+04,
                        4.1732540000e+03, 1.3614770000e+03, 4.9150490000e+02,
                        1.9168570000e+02, 7.9454990000e+01, 3.4549680000e+01,
                        1.5586730000e+01, 7.2354440000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.3837390000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3697280000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.2508900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7481900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1920900000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.4000000000e-05, 5.6200000000e-04, 3.0850000000e-03,
                        1.2298000000e-02, 3.8625000000e-02, 9.9815000000e-02},
              doubles_t{6.7355560000e+02, 1.5958450000e+02, 5.1801700000e+01,
                        1.9673430000e+01, 8.1753450000e+00, 3.6161310000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6734140000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.8924600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7233900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7403800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.8681000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6780800000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1992000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3330000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5541700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3142400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0080730000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6700300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.9584600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8993600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0632000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6930000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.0541700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7792530000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1580220000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2220000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(8): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{6.0000000000e-06, 4.3000000000e-05, 2.2700000000e-04,
                        9.5700000000e-04, 3.4720000000e-03, 1.1196000000e-02,
                        3.2375000000e-02, 8.2835000000e-02, 1.7966800000e-01,
                        3.0520100000e-01, 3.4088900000e-01, 1.7742200000e-01,
                        2.0470000000e-02, -9.5900000000e-04, 1.0530000000e-03,
                        -9.6000000000e-05},
              doubles_t{5.7097980000e+05, 8.5496900000e+04, 1.9456650000e+04,
                        5.5110300000e+03, 1.7979190000e+03, 6.4906770000e+02,
                        2.5313860000e+02, 1.0493790000e+02, 4.5650360000e+01,
                        2.0617620000e+01, 9.5867590000e+00, 4.4925100000e+00,
                        1.8369430000e+00, 8.3480300000e-01, 3.6572200000e-01,
                        1.5703400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4925100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8369430000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3480300000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.0000000000e-06, -1.0000000000e-05, -5.1000000000e-05,
                        -2.1700000000e-04, -7.9100000000e-04, -2.5590000000e-03,
                        -7.5300000000e-03, -1.9782000000e-02, -4.6086000000e-02,
                        -8.9189000000e-02, -1.3755000000e-01, -1.0510800000e-01,
                        1.4478400000e-01, 4.4160000000e-01, 4.4454500000e-01,
                        1.2823200000e-01},
              doubles_t{5.7097980000e+05, 8.5496900000e+04, 1.9456650000e+04,
                        5.5110300000e+03, 1.7979190000e+03, 6.4906770000e+02,
                        2.5313860000e+02, 1.0493790000e+02, 4.5650360000e+01,
                        2.0617620000e+01, 9.5867590000e+00, 4.4925100000e+00,
                        1.8369430000e+00, 8.3480300000e-01, 3.6572200000e-01,
                        1.5703400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6572200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5703400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7326240000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2272000000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.6700000000e-04, 1.4330000000e-03, 7.5480000000e-03,
                        2.8601000000e-02, 8.4361000000e-02, 1.8742700000e-01,
                        2.9795300000e-01, 3.3862300000e-01, 2.4723300000e-01,
                        7.0024000000e-02},
              doubles_t{5.2576380000e+02, 1.2463060000e+02, 4.0343000000e+01,
                        1.5178310000e+01, 6.2450390000e+00, 2.7326240000e+00,
                        1.2272000000e+00, 5.4919000000e-01, 2.4177600000e-01,
                        1.0250500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.4919000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4177600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0250500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7002000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.3625000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4684000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4120000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8000000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.2537000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7084000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1729000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.0800000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.0017000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1158000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.9500000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0207000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5530000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6000000000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
