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
aug-cc-pwcvqz-pp-optri atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(aug_dash_cc_dash_pwcvqz_dash_pp_dash_optri_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(aug_dash_cc_dash_pwcvqz_dash_pp_dash_optri_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("aug-cc-pwcvqz-pp-optri");
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
                                           doubles_t{1.5659560000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.0291410000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1261290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9649900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1883730000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.6068050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2153290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6793800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8901380000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0679320000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9743060000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2590120000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2026860000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4174180000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1755710000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4163280000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8759450000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.6850860000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6787340000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0556140000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.9601210000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.6486840000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3828290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.9956600000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.9191270000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.6819400000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(30): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6806520000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1195790000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.4418920000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1428900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.6314110000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.0447650000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.1993610000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2674230000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9676850000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.7109170000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4171910000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.4823400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.4930720000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3763560000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.7332330000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6111490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8949760000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6267490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.1183980000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2602890000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.3057470000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.6880800000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8620040000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4831900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0353340000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2450240000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(47): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1005430000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.4115940000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.9794570000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1015100000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3122120000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.7574330000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8517490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0859340000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6218240000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0903320000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3316460000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.5557400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2649640000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.5978000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0002450000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.7402400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.7408230000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8239680000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3626490000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.9971000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.3329700000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5578820000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1962280000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3393900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.2121970000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3016100000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(48): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2605080000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.4676120000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.6822680000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0407500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.8327880000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.8919730000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8173450000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1196930000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3539420000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6535200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5143640000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.0107300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.0910250000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.1483240000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.2441100000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.8335400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.8486880000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4547520000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5396980000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.9621200000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.1403530000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.1085640000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3753010000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.8218300000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.6194630000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.5070410000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(79): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.6340260000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0882000000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3692000000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0956420000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1602800000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.7695160000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.2034440000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.2917000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1516340000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.6826050000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.3273000000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.8098400000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0608930000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.0681740000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.8671710000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7735700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9319550000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.1989200000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0522540000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.6795500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9748850000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2035190000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.6163700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.8396500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0514970000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1871650000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(80): {
            shells_t shells;
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.4422510000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.6245240000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3065590000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 0,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{2.9361840000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.2287990000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{8.2223340000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.5047090000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 1,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.8990500000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1838700000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.9081460000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0832480000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.0724700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1388400000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.5841900000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.1482290000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 3,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{6.4992900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.0422490000e+01}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.7428270000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.1655730000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 4,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{5.3331800000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{7.3101260000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.8405280000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.5453900000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 5,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{4.2728700000e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.3166650000e+00}));
            shells.emplace_back(make_shell(pure_t::pure, 6,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{1.3286570000e+00}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
