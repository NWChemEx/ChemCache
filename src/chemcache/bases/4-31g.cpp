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
#include <simde/basis_sets/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
4-31g atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(four_dash_31g_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(four_dash_31g_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("4-31g");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(1): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{3.3494604340e-02, 2.3472695350e-01, 8.1375732610e-01},
              doubles_t{1.8731136960e+01, 2.8253943650e+00, 6.4012169230e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6127775880e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(2): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{4.0139739350e-02, 2.6124609700e-01, 7.9318462460e-01},
              doubles_t{3.8421634000e+01, 5.7780300000e+00, 1.2417740000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9796400000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(5): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7994179600e-02, 1.2469370000e-01,
                                   4.3433537500e-01, 5.6097937400e-01},
                         doubles_t{3.3075285200e+02, 4.9843865000e+01,
                                   1.1117053500e+01, 2.9227243100e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.3038707790e-01, -2.5143438980e-01, 1.2051291990e+00},
              doubles_t{5.6812646210e+00, 1.4544045930e+00, 4.2837857570e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.3742922520e-02, 2.7613305310e-01, 7.7738659620e-01},
              doubles_t{5.6812646210e+00, 1.4544045930e+00, 4.2837857570e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4421917330e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4421917330e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(6): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7725822900e-02, 1.2347786700e-01,
                                   4.3387540000e-01, 5.6150419700e-01},
                         doubles_t{4.8696692800e+02, 7.3371094200e+01,
                                   1.6413457900e+01, 4.3449835600e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.2138374870e-01, -2.2733849750e-01, 1.1851739170e+00},
              doubles_t{8.6735253100e+00, 2.0966192600e+00, 6.0465132900e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.3545384110e-02, 2.9826775710e-01, 7.6210322810e-01},
              doubles_t{8.6735253100e+00, 2.0966192600e+00, 6.0465132900e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8355782980e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8355782980e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(7): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7598251110e-02, 1.2284624110e-01,
                                   4.3378214140e-01, 5.6141821750e-01},
                         doubles_t{6.7127950300e+02, 1.0120166200e+02,
                                   2.2699965900e+01, 6.0406090000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1748929910e-01, -2.1399401620e-01, 1.1745021110e+00},
              doubles_t{1.2393599720e+01, 2.9223828310e+00, 8.3252807680e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.4020344330e-02, 3.1120255510e-01, 7.5274823930e-01},
              doubles_t{1.2393599720e+01, 2.9223828310e+00, 8.3252807680e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2596417300e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2596417300e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(8): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7550627990e-02, 1.2282922300e-01,
                                   4.3488358380e-01, 5.6001080380e-01},
                         doubles_t{8.8327286000e+02, 1.3312928000e+02,
                                   2.9906407900e+01, 7.9786771600e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1340100290e-01, -1.7728646590e-01, 1.1504079290e+00},
              doubles_t{1.6194446640e+01, 3.7800860220e+00, 1.0709835750e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.8545274710e-02, 3.3122543500e-01, 7.3460787810e-01},
              doubles_t{1.6194446640e+01, 3.7800860220e+00, 1.0709835750e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8387984070e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8387984070e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(9): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7475760900e-02, 1.2252308900e-01,
                                   4.3499850200e-01, 5.5981216700e-01},
                         doubles_t{1.1261626900e+03, 1.6974315700e+02,
                                   3.8181511200e+01, 1.0212035900e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1105707950e-01, -1.6832210180e-01, 1.1436255550e+00},
              doubles_t{2.1495366700e+01, 4.9897775700e+00, 1.4035738600e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.9888750800e-02, 3.3938751000e-01, 7.2795898100e-01},
              doubles_t{2.1495366700e+01, 4.9897775700e+00, 1.4035738600e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7303183500e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7303183500e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(10): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.7423805400e-02, 1.2227274500e-01,
                                   4.3501423200e-01, 5.5971464200e-01},
                         doubles_t{1.3979320800e+03, 2.1076978100e+02,
                                   4.7467256900e+01, 1.2722626300e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.0960943870e-01, -1.6412488950e-01, 1.1401515860e+00},
              doubles_t{2.7213033200e+01, 6.2941343500e+00, 1.7600512500e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{7.0440306680e-02, 3.4399304690e-01, 7.2451495980e-01},
              doubles_t{2.7213033200e+01, 6.2941343500e+00, 1.7600512500e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.6186699200e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.6186699200e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(15): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.8521313720e-02, 1.2990486420e-01,
                                   4.5510028860e-01, 5.3313186170e-01},
                         doubles_t{3.0186717800e+03, 4.5512712100e+02,
                                   1.0231473000e+02, 2.7617847300e+01}));
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{-2.4750296130e-02, -1.3509246010e-01,
                                   2.2773608020e-01, 8.7559311690e-01},
                         doubles_t{1.1442940100e+02, 2.6582295900e+01,
                                   7.8718889000e+00, 2.4878572500e+00}));
            shells.emplace_back(
              make_shell(pure_t::pure, 1,
                         doubles_t{2.7414002550e-02, 1.6907914230e-01,
                                   4.6910208990e-01, 5.1815306000e-01},
                         doubles_t{1.1442940100e+02, 2.6582295900e+01,
                                   7.8718889000e+00, 2.4878572500e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-4.5119223380e-02, -8.5047299710e-01, 1.5962858630e+00},
              doubles_t{5.0750619000e+01, 1.6728624200e+00, 6.2109741200e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{3.7790712150e-03, -4.6343840930e-02, 1.0339443000e+00},
              doubles_t{5.0750619000e+01, 1.6728624200e+00, 6.2109741200e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6701600700e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6701600700e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(16): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.8492123620e-02, 1.2982202210e-01,
                                   4.5504178740e-01, 5.3300835650e-01},
                         doubles_t{3.4421244100e+03, 5.1891310000e+02,
                                   1.1669090300e+02, 3.1571647200e+01}));
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{-2.7264610620e-02, -1.4248341510e-01,
                                   2.5970435220e-01, 8.5254729550e-01},
                         doubles_t{1.2744057600e+02, 2.9747667300e+01,
                                   8.8346642800e+00, 2.8173898200e+00}));
            shells.emplace_back(
              make_shell(pure_t::pure, 1,
                         doubles_t{2.9151999540e-02, 1.7795967630e-01,
                                   4.8362371270e-01, 4.9425530270e-01},
                         doubles_t{1.2744057600e+02, 2.9747667300e+01,
                                   8.8346642800e+00, 2.8173898200e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.7753152500e-01, -4.5764345830e-01, 1.4316842800e+00},
              doubles_t{3.7291853700e+00, 1.4067701700e+00, 5.4810996900e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-3.3750926340e-02, 1.4571104520e-01, 8.9828874420e-01},
              doubles_t{3.7291853700e+00, 1.4067701700e+00, 5.4810996900e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7038090500e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7038090500e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(17): {
            shells_t shells;
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{1.8379431120e-02, 1.2914012320e-01,
                                   4.5404489060e-01, 5.3443943670e-01},
                         doubles_t{3.9103026900e+03, 5.8955180700e+02,
                                   1.3259392400e+02, 3.5903542500e+01}));
            shells.emplace_back(
              make_shell(pure_t::pure, 0,
                         doubles_t{-2.6743323030e-02, -1.4469118220e-01,
                                   2.5170356930e-01, 8.5982038190e-01},
                         doubles_t{1.4776535300e+02, 3.4506075300e+01,
                                   1.0286471500e+01, 3.3111473800e+00}));
            shells.emplace_back(
              make_shell(pure_t::pure, 1,
                         doubles_t{2.8864468810e-02, 1.7796467010e-01,
                                   4.8699980720e-01, 4.8901845020e-01},
                         doubles_t{1.4776535300e+02, 3.4506075300e+01,
                                   1.0286471500e+01, 3.3111473800e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.7039627540e-01, -3.4162971950e-01, 1.3500244820e+00},
              doubles_t{4.2802849100e+00, 1.6410166700e+00, 6.1447850300e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-3.6702885140e-02, 1.9184924220e-01, 8.6433768180e-01},
              doubles_t{4.2802849100e+00, 1.6410166700e+00, 6.1447850300e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9565941100e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9565941100e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
