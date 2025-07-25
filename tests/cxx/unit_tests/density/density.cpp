/*
 * Copyright 2025 NWChemEx-Project
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
#include <simde/chemical_system/atomic_density_matrix.hpp>
#include <simde/density/initial_density.hpp>
#include <simde/types.hpp>

using atomic_dm_pt        = simde::AtomicDensityMatrixFromZ;
using sad_rho_pt          = simde::InitialDensity;
using molecular_basis_pt  = simde::MolecularBasisSet;
using atomic_number_t     = simde::type::atomic_number;
using molecule_t          = simde::type::molecule;
using atom_t              = simde::type::atom;
using hamiltonian_t       = simde::type::hamiltonian;
using V_en_t              = simde::type::V_en_type;
using aos_t               = simde::type::aos;
using udouble_t           = tensorwrapper::types::udouble;
using shape_t             = tensorwrapper::shape::Smooth;
using double_allocator_t  = tensorwrapper::allocator::Eigen<double>;
using double_il_t         = typename double_allocator_t::rank2_il;
using udouble_allocator_t = tensorwrapper::allocator::Eigen<udouble_t>;
using udouble_il_t        = typename udouble_allocator_t::rank2_il;

using Catch::Matchers::Message;

TEST_CASE("sto-3g atomic density matrix") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_dm_mod = mm.at("sto-3g atomic density matrix");
    auto rt           = mm.get_runtime();

    SECTION("Hydrogen") {
        atomic_number_t Z = 1;
        shape_t shape{1, 1};
        SECTION("Exact Values") {
            const auto& rv = atom_dm_mod.run_as<atomic_dm_pt>(Z);
            double_allocator_t alloc(rt);
            auto buffer = alloc.construct(double_il_t{{1.0}});
            simde::type::tensor corr(shape, std::move(buffer));
            REQUIRE(rv == corr);
        }
        SECTION("Uncertain Values") {
            atom_dm_mod.change_input("With UQ?", true);
            const auto& rv = atom_dm_mod.run_as<atomic_dm_pt>(Z);
            udouble_allocator_t alloc(rt);
            auto buffer = alloc.construct(udouble_il_t{{1.0}});
            simde::type::tensor corr(shape, std::move(buffer));
            REQUIRE(rv == corr);
        }
    }

    SECTION("Out of Range") {
        atomic_number_t Z = 1000;
        REQUIRE_THROWS_MATCHES(atom_dm_mod.run_as<atomic_dm_pt>(Z),
                               std::out_of_range,
                               Message("Atomic Density not available for Z"));
    }
}

TEST_CASE("SAD Density") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& sad_mod = mm.at("sto-3g SAD density");
    auto rt       = mm.get_runtime();

    atom_t h1{"H", 1, 0.0, 1.0, 1.0, 1.0};
    atom_t h2{"H", 1, 0.0, 2.0, 2.0, 2.0};
    molecule_t mol{h1, h2};
    hamiltonian_t H;
    H.emplace_back(
      1.0, std::make_unique<V_en_t>(mol.electrons(), mol.nuclei().as_nuclei()));

    // Correct values: Density Matrix (2x2 identity) and AOs
    shape_t shape{2, 2};
    double_allocator_t alloc(rt);
    auto buffer = alloc.construct(double_il_t{{0.5, 0.0}, {0.0, 0.5}});
    simde::type::tensor corr_value(shape, std::move(buffer));

    auto& mol_bs_mod = mm.at("sto-3g");
    auto aobs        = mol_bs_mod.run_as<molecular_basis_pt>(mol);
    aos_t corr_aos(aobs);

    // Run the module
    const auto& rv = sad_mod.run_as<sad_rho_pt>(H);
    REQUIRE(rv.value() == corr_value);
    REQUIRE(rv.basis_set() == corr_aos);
}