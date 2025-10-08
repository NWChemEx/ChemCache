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

#include "../catch.hpp"
#include "chemcache/chemcache.hpp"
#include <simde/basis_set/basis_set.hpp>
#include <simde/types.hpp>

using atomic_basis_pt    = simde::AtomicBasisSetFromZ;
using atomic_basis_t     = simde::type::atomic_basis_set;
using cg_t               = simde::type::contracted_gaussian;
using molecular_basis_pt = simde::MolecularBasisSet;
using molecular_basis_t  = simde::type::ao_basis_set;
using molecule_t         = simde::type::molecule;
using atom_t             = simde::type::atom;
using atomic_number_t    = simde::type::atomic_number;
using doubles_t          = std::vector<double>;

using Catch::Matchers::Message;

TEST_CASE("Atomic Basis Set") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_bs_mod = mm.at("sto-3g atomic basis");

    SECTION("sto-3g Hydrogen") {
        atomic_number_t Z = 1;
        auto rv           = atom_bs_mod.run_as<atomic_basis_pt>(Z);
        doubles_t cs{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01};
        doubles_t es{3.4252509140e+00, 6.2391372980e-01, 1.6885540400e-01};
        cg_t cg(cs.begin(), cs.end(), es.begin(), es.end(), 0.0, 0.0, 0.0);

        atomic_basis_t corr("sto-3g", 1, 0.0, 0.0, 0.0);
        corr.add_shell(chemist::ShellType::pure, 0, cg);
        REQUIRE(rv == corr);
    }

    SECTION("sto-3g Oxygen") {
        atomic_number_t Z = 8;
        auto rv           = atom_bs_mod.run_as<atomic_basis_pt>(Z);
        doubles_t cs0{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01};
        doubles_t es0{1.3070932140e+02, 2.3808866050e+01, 6.4436083130e+00};
        cg_t cg0(cs0.begin(), cs0.end(), es0.begin(), es0.end(), 0.0, 0.0, 0.0);
        doubles_t cs1{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01};
        doubles_t es1{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01};
        cg_t cg1(cs1.begin(), cs1.end(), es1.begin(), es1.end(), 0.0, 0.0, 0.0);
        doubles_t cs2{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01};
        doubles_t es2{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01};
        cg_t cg2(cs2.begin(), cs2.end(), es2.begin(), es2.end(), 0.0, 0.0, 0.0);

        atomic_basis_t corr("sto-3g", 8, 0.0, 0.0, 0.0);
        corr.add_shell(chemist::ShellType::pure, 0, cg0);
        corr.add_shell(chemist::ShellType::pure, 0, cg1);
        corr.add_shell(chemist::ShellType::pure, 1, cg2);
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        atomic_number_t Z = 1000;
        auto z_string     = std::to_string(Z);
        REQUIRE_THROWS_MATCHES(
          atom_bs_mod.run_as<atomic_basis_pt>(Z), std::out_of_range,
          Message("Basis Set not available for Z: " + z_string));
    }
}

TEST_CASE("Molecular Basis Set") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& mol_bs_mod = mm.at("sto-3g");

    atom_t h1{"H", 1, 0.0, 1.0, 1.0, 1.0};
    atom_t h2{"H", 1, 0.0, 2.0, 2.0, 2.0};
    molecule_t in{h1, h2};

    auto h_aos = [](auto r) {
        doubles_t cs{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01};
        doubles_t es{3.4252509140e+00, 6.2391372980e-01, 1.6885540400e-01};
        cg_t cg(cs.begin(), cs.end(), es.begin(), es.end(), 0.0, 0.0, 0.0);

        atomic_basis_t aos("sto-3g", 1, r, r, r);
        aos.add_shell(chemist::ShellType::pure, 0, cg);
        return aos;
    };

    auto atom_mod = pluginplay::make_lambda<atomic_basis_pt>([=](auto&& Z) {
        REQUIRE(Z == 1);
        return h_aos(0.0);
    });
    mol_bs_mod.change_submod("Atomic Basis", atom_mod);

    molecular_basis_t corr;
    corr.add_center(h_aos(1.0));
    corr.add_center(h_aos(2.0));

    auto rv = mol_bs_mod.run_as<molecular_basis_pt>(in);
    REQUIRE(rv == corr);
}
