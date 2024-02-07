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
sto-3g* atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(sto_dash_3g_star_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(sto_dash_3g_star_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    // Basis Set name and origin point
    std::string name("sto-3g*");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    switch(Z) {
        case(11): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{2.5077243000e+02, 4.5678511170e+01, 1.2362387760e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{1.2040192740e+01, 2.7978818590e+00, 9.0995801700e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{1.2040192740e+01, 2.7978818590e+00, 9.0995801700e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{1.4787406220e+00, 4.1256488010e-01, 1.6147509790e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{1.4787406220e+00, 4.1256488010e-01, 1.6147509790e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.0000000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(12): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{2.9923741370e+02, 5.4506468450e+01, 1.4751577520e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{1.5121823520e+01, 3.5139865790e+00, 1.1428574980e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{1.5121823520e+01, 3.5139865790e+00, 1.1428574980e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{9.0000000000e-02}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(13): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{3.5142147670e+02, 6.4011860670e+01, 1.7324107610e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{1.8899396210e+01, 4.3918132330e+00, 1.4283539700e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{1.8899396210e+01, 4.3918132330e+00, 1.4283539700e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{1.3954482930e+00, 3.8932653180e-01, 1.5237976590e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(14): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{4.0779755140e+02, 7.4280833050e+01, 2.0103292290e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{2.3193656060e+01, 5.3897068710e+00, 1.7528999520e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{2.3193656060e+01, 5.3897068710e+00, 1.7528999520e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{1.4787406220e+00, 4.1256488010e-01, 1.6147509790e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{1.4787406220e+00, 4.1256488010e-01, 1.6147509790e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(15): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{4.6836563780e+02, 8.5313385590e+01, 2.3089131560e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{2.8032639580e+01, 6.5141825770e+00, 2.1186143520e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{2.8032639580e+01, 6.5141825770e+00, 2.1186143520e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{1.7431032310e+00, 4.8632137710e-01, 1.9034289090e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{1.7431032310e+00, 4.8632137710e-01, 1.9034289090e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(16): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{5.3312573590e+02, 9.7109518300e+01, 2.6281625420e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{3.3329751730e+01, 7.7451175210e+00, 2.5189525990e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{3.3329751730e+01, 7.7451175210e+00, 2.5189525990e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{2.0291942740e+00, 5.6614005180e-01, 2.2158337920e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{2.0291942740e+00, 5.6614005180e-01, 2.2158337920e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(17): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{6.0134561360e+02, 1.0953585420e+02, 2.9644676860e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{3.8960418890e+01, 9.0535634770e+00, 2.9444998340e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{3.8960418890e+01, 9.0535634770e+00, 2.9444998340e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{2.1293864950e+00, 5.9409342740e-01, 2.3252414100e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{2.1293864950e+00, 5.9409342740e-01, 2.3252414100e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        case(18): {
            shells_t shells;
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01},
              doubles_t{6.7444651840e+02, 1.2285127530e+02, 3.3248349450e+01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01},
              doubles_t{4.5164243920e+01, 1.0495199000e+01, 3.4133644480e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01},
              doubles_t{4.5164243920e+01, 1.0495199000e+01, 3.4133644480e+00}));
            shells.emplace_back(make_shell(
              pure_t::pure, 0,
              doubles_t{-2.1962036900e-01, 2.2559543360e-01, 9.0039842600e-01},
              doubles_t{2.6213665180e+00, 7.3135460500e-01, 2.8624723560e-01}));
            shells.emplace_back(make_shell(
              pure_t::pure, 1,
              doubles_t{1.0587604290e-02, 5.9516700530e-01, 4.6200101200e-01},
              doubles_t{2.6213665180e+00, 7.3135460500e-01, 2.8624723560e-01}));
            shells.emplace_back(make_shell(pure_t::pure, 2,
                                           doubles_t{1.0000000000e+00},
                                           doubles_t{3.9000000000e-01}));
            abs_t atom_bs(name, Z, r0, shells.begin(), shells.end());
            return atomic_basis_pt::wrap_results(rv, atom_bs);
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
