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
#include <simde/basis_set/atomic_basis_set.hpp>
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
6-31+g atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(six_dash_31_plus_g_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(six_dash_31_plus_g_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("6-31+g");
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
        case(3): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{2.1426078100e-03, 1.6208871500e-02, 7.7315572500e-02,
                        2.4578605200e-01, 4.7018900400e-01, 3.4547084500e-01},
              doubles_t{6.4241891500e+02, 9.6798515300e+01, 2.2091121200e+01,
                        6.2010702500e+00, 1.9351176800e+00, 6.3673578900e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-3.5091745740e-02, -1.9123284310e-01, 1.0839877950e+00},
              doubles_t{2.3249184080e+00, 6.3243035560e-01, 7.9053434750e-02}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{8.9415080430e-03, 1.4100946400e-01, 9.4536369530e-01},
              doubles_t{2.3249184080e+00, 6.3243035560e-01, 7.9053434750e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.5961971750e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.5961971750e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4000000000e-03}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4000000000e-03}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(4): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.9447575900e-03, 1.4835052000e-02, 7.2090546290e-02,
                        2.3715415000e-01, 4.6919865190e-01, 3.5652022790e-01},
              doubles_t{1.2645856900e+03, 1.8993680600e+02, 4.3159089000e+01,
                        1.2098662700e+01, 3.8063232200e+00, 1.2728903000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1264872850e-01, -2.2950640790e-01, 1.1869167640e+00},
              doubles_t{3.1964630980e+00, 7.4781330380e-01, 2.1996633020e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{5.5980199800e-02, 2.6155061100e-01, 7.9397233890e-01},
              doubles_t{3.1964630980e+00, 7.4781330380e-01, 2.1996633020e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.2309900700e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.2309900700e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0700000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0700000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(5): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8662745900e-03, 1.4251481700e-02, 6.9551618500e-02,
                        2.3257293300e-01, 4.6707871200e-01, 3.6343144000e-01},
              doubles_t{2.0688822500e+03, 3.1064957000e+02, 7.0683033000e+01,
                        1.9861080300e+01, 6.2993048400e+00, 2.1270269700e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.3039379740e-01, -1.3078895140e-01, 1.1309444840e+00},
              doubles_t{4.7279710710e+00, 1.1903377360e+00, 3.5941168290e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{7.4597579920e-02, 3.0784667710e-01, 7.4345683420e-01},
              doubles_t{4.7279710710e+00, 1.1903377360e+00, 3.5941168290e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2675124690e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2675124690e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1500000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1500000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(6): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8347371320e-03, 1.4037322810e-02, 6.8842622260e-02,
                        2.3218444320e-01, 4.6794134840e-01, 3.6231198530e-01},
              doubles_t{3.0475248800e+03, 4.5736951800e+02, 1.0394868500e+02,
                        2.9210155300e+01, 9.2866629600e+00, 3.1639269600e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1933241980e-01, -1.6085415170e-01, 1.1434564380e+00},
              doubles_t{7.8682723500e+00, 1.8812885400e+00, 5.4424925800e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.8999066590e-02, 3.1642396100e-01, 7.4430829090e-01},
              doubles_t{7.8682723500e+00, 1.8812885400e+00, 5.4424925800e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6871447820e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6871447820e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3800000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3800000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(7): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8347721600e-03, 1.3994627000e-02, 6.8586551810e-02,
                        2.3224087300e-01, 4.6906994810e-01, 3.6045519910e-01},
              doubles_t{4.1735114600e+03, 6.2745791100e+02, 1.4290209300e+02,
                        4.0234329300e+01, 1.2820212900e+01, 4.3904370100e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1496118170e-01, -1.6911747860e-01, 1.1458519470e+00},
              doubles_t{1.1626361860e+01, 2.7162798070e+00, 7.7221839660e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{6.7579743880e-02, 3.2390729590e-01, 7.4089513980e-01},
              doubles_t{1.1626361860e+01, 2.7162798070e+00, 7.7221839660e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1203149750e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1203149750e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3900000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3900000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(8): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8310744300e-03, 1.3950172200e-02, 6.8445078100e-02,
                        2.3271433600e-01, 4.7019289800e-01, 3.5852085300e-01},
              doubles_t{5.4846716600e+03, 8.2523494600e+02, 1.8804695800e+02,
                        5.2964500000e+01, 1.6897570400e+01, 5.7996353400e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.1077754950e-01, -1.4802626270e-01, 1.1307670150e+00},
              doubles_t{1.5539616250e+01, 3.5999335860e+00, 1.0137617500e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{7.0874268230e-02, 3.3975283910e-01, 7.2715857730e-01},
              doubles_t{1.5539616250e+01, 3.5999335860e+00, 1.0137617500e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7000582260e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7000582260e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4500000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4500000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(9): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8196169010e-03, 1.3916079610e-02, 6.8405324530e-02,
                        2.3318576010e-01, 4.7126743920e-01, 3.5661854620e-01},
              doubles_t{7.0017130900e+03, 1.0513660900e+03, 2.3928569000e+02,
                        6.7397445300e+01, 2.1519957300e+01, 7.4031013000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.0850697510e-01, -1.4645165810e-01, 1.1286885810e+00},
              doubles_t{2.0847952800e+01, 4.8083083400e+00, 1.3440698600e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{7.1628724240e-02, 3.4591210270e-01, 7.2246995640e-01},
              doubles_t{2.0847952800e+01, 4.8083083400e+00, 1.3440698600e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.5815139300e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.5815139300e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0760000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0760000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(10): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8843480500e-03, 1.4336899400e-02, 7.0109623310e-02,
                        2.3737326600e-01, 4.7300712610e-01, 3.4840124100e-01},
              doubles_t{8.4258515300e+03, 1.2685194000e+03, 2.8962141400e+02,
                        8.1859004000e+01, 2.6251507900e+01, 9.0947205100e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-1.0711828720e-01, -1.4616382130e-01, 1.1277735030e+00},
              doubles_t{2.6532131000e+01, 6.1017550100e+00, 1.6962715300e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{7.1909588510e-02, 3.4951337200e-01, 7.1994051210e-01},
              doubles_t{2.6532131000e+01, 6.1017550100e+00, 1.6962715300e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4581870000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4581870000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3000000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(11): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.9376592770e-03, 1.4806994480e-02, 7.2705472880e-02,
                        2.5262890580e-01, 4.9324181600e-01, 3.1316888320e-01},
              doubles_t{9.9932000000e+03, 1.4998900000e+03, 3.4195100000e+02,
                        9.4679600000e+01, 2.9734500000e+01, 1.0006300000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-3.5420835040e-03, -4.3958843480e-02, -1.0975210860e-01,
                        1.8739818540e-01, 6.4669963970e-01, 3.0605830270e-01},
              doubles_t{1.5096300000e+02, 3.5587800000e+01, 1.1168300000e+01,
                        3.9020100000e+00, 1.3817700000e+00, 4.6638200000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{5.0016597100e-03, 3.5510897940e-02, 1.4282499170e-01,
                        3.3861998030e-01, 4.5157897380e-01, 2.7327098410e-01},
              doubles_t{1.5096300000e+02, 3.5587800000e+01, 1.1168300000e+01,
                        3.9020100000e+00, 1.3817700000e+00, 4.6638200000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.4850315930e-01, -1.3170408440e-01, 1.2335207910e+00},
              doubles_t{4.9796600000e-01, 8.4352900000e-02, 6.6635000000e-02}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-2.3022500430e-02, 9.5035901760e-01, 5.9857901110e-02},
              doubles_t{4.9796600000e-01, 8.4352900000e-02, 6.6635000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5954400000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5954400000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6000000000e-03}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6000000000e-03}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(12): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.9778293170e-03, 1.5113994780e-02, 7.3910774480e-02,
                        2.4919091400e-01, 4.8792783160e-01, 3.1966188960e-01},
              doubles_t{1.1722800000e+04, 1.7599300000e+03, 4.0084600000e+02,
                        1.1280700000e+02, 3.5999700000e+01, 1.2182800000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-3.2371704710e-03, -4.1007905970e-02, -1.1260001640e-01,
                        1.4863302160e-01, 6.1649708980e-01, 3.6482905310e-01},
              doubles_t{1.8918000000e+02, 4.5211900000e+01, 1.4356300000e+01,
                        5.1388600000e+00, 1.9065200000e+00, 7.0588700000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.9281299210e-03, 3.4988799440e-02, 1.4072499770e-01,
                        3.3364199470e-01, 4.4493999290e-01, 2.6925399570e-01},
              doubles_t{1.8918000000e+02, 4.5211900000e+01, 1.4356300000e+01,
                        5.1388600000e+00, 1.9065200000e+00, 7.0588700000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1229089850e-01, -1.0798545700e-01, 1.1758449770e+00},
              doubles_t{9.2934000000e-01, 2.6903500000e-01, 1.1737900000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-2.2419181230e-02, 1.9227083900e-01, 8.4618029160e-01},
              doubles_t{9.2934000000e-01, 2.6903500000e-01, 1.1737900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2106100000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2106100000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4600000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4600000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(13): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.9426699470e-03, 1.4859899590e-02, 7.2849398000e-02,
                        2.4682999320e-01, 4.8725798660e-01, 3.2349599110e-01},
              doubles_t{1.3983100000e+04, 2.0987500000e+03, 4.7770500000e+02,
                        1.3436000000e+02, 4.2870900000e+01, 1.4518900000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.9261900280e-03, -3.7408300360e-02, -1.1448700110e-01,
                        1.1563500110e-01, 6.1259500580e-01, 3.9379900370e-01},
              doubles_t{2.3966800000e+02, 5.7441900000e+01, 1.8285900000e+01,
                        6.5991400000e+00, 2.4904900000e+00, 9.4454500000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.6028455820e-03, 3.3198968130e-02, 1.3628186920e-01,
                        3.3047568280e-01, 4.4914556890e-01, 2.6570374500e-01},
              doubles_t{2.3966800000e+02, 5.7441900000e+01, 1.8285900000e+01,
                        6.5991400000e+00, 2.4904900000e+00, 9.4454500000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.2760692450e-01, 1.4458358730e-03, 1.0927944390e+00},
              doubles_t{1.2779000000e+00, 3.9759000000e-01, 1.6009500000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.7512601890e-02, 2.4453302640e-01, 8.0493408670e-01},
              doubles_t{1.2779000000e+00, 3.9759000000e-01, 1.6009500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5657700000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5657700000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1800000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1800000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(14): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.9594802160e-03, 1.4928801640e-02, 7.2847808010e-02,
                        2.4613002710e-01, 4.8591405350e-01, 3.2500203580e-01},
              doubles_t{1.6115900000e+04, 2.4255800000e+03, 5.5386700000e+02,
                        1.5634000000e+02, 5.0068300000e+01, 1.7017800000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.7809414150e-03, -3.5714618170e-02, -1.1498505850e-01,
                        9.3563447600e-02, 6.0301730680e-01, 4.1895921310e-01},
              doubles_t{2.9271800000e+02, 6.9873100000e+01, 2.2336300000e+01,
                        8.1503900000e+00, 3.1345800000e+00, 1.2254300000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.4382645210e-03, 3.2667933280e-02, 1.3472113720e-01,
                        3.2867833480e-01, 4.4964045800e-01, 2.6137226620e-01},
              doubles_t{2.9271800000e+02, 6.9873100000e+01, 2.2336300000e+01,
                        8.1503900000e+00, 3.1345800000e+00, 1.2254300000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.4463100420e-01, 4.3157377170e-03, 1.0981845080e+00},
              doubles_t{1.7273800000e+00, 5.7292200000e-01, 2.2219200000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.7795106050e-02, 2.5353908630e-01, 8.0066927240e-01},
              doubles_t{1.7273800000e+00, 5.7292200000e-01, 2.2219200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7836900000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7836900000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.3100000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.3100000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(15): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8515989230e-03, 1.4206191740e-02, 6.9999459280e-02,
                        2.4007886030e-01, 4.8476171800e-01, 3.3519980500e-01},
              doubles_t{1.9413300000e+04, 2.9094200000e+03, 6.6136400000e+02,
                        1.8575900000e+02, 5.9194300000e+01, 2.0031000000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.7821701050e-03, -3.6049901350e-02, -1.1663100440e-01,
                        9.6832803640e-02, 6.1441802310e-01, 4.0379801520e-01},
              doubles_t{3.3947800000e+02, 8.1010100000e+01, 2.5878000000e+01,
                        9.4522100000e+00, 3.6656600000e+00, 1.4674600000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.5646161910e-03, 3.3693571880e-02, 1.3975488340e-01,
                        3.3936171680e-01, 4.5092062370e-01, 2.3858580090e-01},
              doubles_t{3.3947800000e+02, 8.1010100000e+01, 2.5878000000e+01,
                        9.4522100000e+00, 3.6656600000e+00, 1.4674600000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.5292411390e-01, 3.2851844680e-02, 1.0812547620e+00},
              doubles_t{2.1562300000e+00, 7.4899700000e-01, 2.8314500000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.7765312730e-02, 2.7405819640e-01, 7.8542156300e-01},
              doubles_t{2.1562300000e+00, 7.4899700000e-01, 2.8314500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.9831700000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.9831700000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4800000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4800000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(16): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8692408490e-03, 1.4230306460e-02, 6.9696231660e-02,
                        2.3848710830e-01, 4.8330721950e-01, 3.3807415360e-01},
              doubles_t{2.1917100000e+04, 3.3014900000e+03, 7.5414600000e+02,
                        2.1271100000e+02, 6.7989600000e+01, 2.3051500000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.3767704990e-03, -3.1693006650e-02, -1.1331702380e-01,
                        5.6090011770e-02, 5.9225512430e-01, 4.5500609550e-01},
              doubles_t{4.2373500000e+02, 1.0071000000e+02, 3.2159900000e+01,
                        1.1807900000e+01, 4.6311000000e+00, 1.8702500000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{4.0610099820e-03, 3.0681299860e-02, 1.3045199940e-01,
                        3.2720499850e-01, 4.5285099800e-01, 2.5604199890e-01},
              doubles_t{4.2373500000e+02, 1.0071000000e+02, 3.2159900000e+01,
                        1.1807900000e+01, 4.6311000000e+00, 1.8702500000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.5037311420e-01, 6.6956763100e-02, 1.0545062690e+00},
              doubles_t{2.6158400000e+00, 9.2216700000e-01, 3.4128700000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.4510489550e-02, 3.1026277650e-01, 7.5448245650e-01},
              doubles_t{2.6158400000e+00, 9.2216700000e-01, 3.4128700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1716700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1716700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0500000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0500000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(17): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8329598480e-03, 1.4034198830e-02, 6.9097394260e-02,
                        2.3745198030e-01, 4.8303395990e-01, 3.3985597180e-01},
              doubles_t{2.5180100000e+04, 3.7803500000e+03, 8.6047400000e+02,
                        2.4214500000e+02, 7.7334900000e+01, 2.6247000000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.2973914170e-03, -3.0713718940e-02, -1.1252806940e-01,
                        4.5016327760e-02, 5.8935336340e-01, 4.6520628680e-01},
              doubles_t{4.9176500000e+02, 1.1698400000e+02, 3.7415300000e+01,
                        1.3783400000e+01, 5.4521500000e+00, 2.2258800000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{3.9894008790e-03, 3.0317706680e-02, 1.2988002860e-01,
                        3.2795107230e-01, 4.5352710000e-01, 2.5215405560e-01},
              doubles_t{4.9176500000e+02, 1.1698400000e+02, 3.7415300000e+01,
                        1.3783400000e+01, 5.4521500000e+00, 2.2258800000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.5182802800e-01, 6.1589251410e-02, 1.0601843280e+00},
              doubles_t{3.1864900000e+00, 1.1442700000e+00, 4.2037700000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.4299314720e-02, 3.2357233310e-01, 7.4350776530e-01},
              doubles_t{3.1864900000e+00, 1.1442700000e+00, 4.2037700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4265700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4265700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8300000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8300000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(18): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.8252601920e-03, 1.3968601470e-02, 6.8707307230e-02,
                        2.3620402490e-01, 4.8221405080e-01, 3.4204303600e-01},
              doubles_t{2.8348300000e+04, 4.2576200000e+03, 9.6985700000e+02,
                        2.7326300000e+02, 8.7369500000e+01, 2.9686700000e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1597208950e-03, -2.9077512060e-02, -1.1082704600e-01,
                        2.7699911480e-02, 5.7761323950e-01, 4.8868820260e-01},
              doubles_t{5.7589100000e+02, 1.3681600000e+02, 4.3809800000e+01,
                        1.6209400000e+01, 6.4608400000e+00, 2.6511400000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{3.8066498420e-03, 2.9230498790e-02, 1.2646699480e-01,
                        3.2350998660e-01, 4.5489598110e-01, 2.5662998940e-01},
              doubles_t{5.7589100000e+02, 1.3681600000e+02, 4.3809800000e+01,
                        1.6209400000e+01, 6.4608400000e+00, 2.6511400000e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.5559296040e-01, 3.7806742060e-02, 1.0805640600e+00},
              doubles_t{3.8602800000e+00, 1.4137300000e+00, 5.1664600000e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{-1.5919690400e-02, 3.2464580420e-01, 7.4398955120e-01},
              doubles_t{3.8602800000e+00, 1.4137300000e+00, 5.1664600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7388800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7388800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.0000000000e-02}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.0000000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
