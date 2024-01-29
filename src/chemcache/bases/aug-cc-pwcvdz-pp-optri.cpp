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
aug-cc-pwcvdz-pp-optri atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(aug_dash_cc_dash_pwcvdz_dash_pp_dash_optri_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(aug_dash_cc_dash_pwcvdz_dash_pp_dash_optri_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("aug-cc-pwcvdz-pp-optri");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(29): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3615150000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4436000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3742140000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4483190000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4686000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1092380000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.3506280000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0277600000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.8788100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.3122300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5688530000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4603340000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1723270000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0688050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.1144700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1890890000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5693780000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3523440000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9762600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4168310000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3687560000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0874870000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5582100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.0312760000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6281350000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.4472130000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(30): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3700190000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.1444660000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2237420000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4850770000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1816000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3262550000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5673670000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9247800000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4626400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3571500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8613380000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5009300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9823310000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6151700000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0748830000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5132730000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3711850000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6204940000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.5101300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7163670000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4469590000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6776790000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0182460000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.1262860000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0970450000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4053200000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(47): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2639210000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.5539920000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6093170000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.9113700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9950200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1469780000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7147390000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.1799100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7497360000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.2696900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1001460000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3441670000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2406120000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3016400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5366500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1659440000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.3985980000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6532840000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.6316400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.0187400000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0883250000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3945630000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8335200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1015640000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3194290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2667700000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(48): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2869320000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.3258540000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5123050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.9460800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4276000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2571740000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.5174520000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7467910000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8560780000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.2131100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1832180000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.3798510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9417370000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0469140000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.5284000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5377210000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.8405300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9481920000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7311600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.1320870000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4520730000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5944050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.8856900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4499030000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5154460000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.5575350000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(79): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5054590000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0041750000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.7854040000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2661230000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4507100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.7283290000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.0880680000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.0522510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5092490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7702100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0167350000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.8013230000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.5419460000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3564320000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.0969300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.3994110000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.6716510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1253190000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.6652600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.4745740000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9860180000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3999910000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7447100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9044630000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2824050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7224630000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(80): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3373940000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.9239300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7967510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4314090000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2733300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2115300000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1483160000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6630510000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4580420000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4860000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0679550000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.1339590000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.7504670000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4854680000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9450900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3062810000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5539610000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7121290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.5591400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8929370000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2679300000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5760520000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.6214800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9810210000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2624060000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8808450000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
