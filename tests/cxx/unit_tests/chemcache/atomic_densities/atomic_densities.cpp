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

#include "chemcache/chemcache.hpp"
#include <catch2/catch.hpp>
#include <simde/simde.hpp>

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using atomic_basis_t  = simde::type::atomic_basis_set;
using atomic_den_pt   = simde::AtomDenFromZ;
using atomic_den_t    = simde::type::el_density;
using tensor_t        = simde::type::tensor;
using ao_bs_t         = simde::type::ao_basis_set;
using aos_t           = simde::type::ao_space;

using Catch::Matchers::Message;

TEST_CASE("Atomic Densities") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_dm_mod = mm.at("sto-3g atomic dm");

    SECTION("sto-3g Hydrogen") {
        auto rho = atom_dm_mod.run_as<atomic_den_pt>(1ul);
        tensor_t corr_dm{1.00000000};
        atomic_basis_t corr_bs("sto-3g", 1, 0.0, 0.0, 0.0);
        corr_bs.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{1.5432896730e-01, 5.3532814230e-01,
                                           4.4463454220e-01},
                       std::vector<double>{3.4252509140e+00, 6.2391372980e-01,
                                           1.6885540400e-01});
        ao_bs_t ao_bs;
        ao_bs.add_center(corr_bs);
        aos_t corr_aos(ao_bs);

        atomic_den_t corr_rho(corr_dm, corr_aos);

        REQUIRE(rho == corr_rho);
    }

    SECTION("sto-3g Oxygen") {
        auto rho = atom_dm_mod.run_as<atomic_den_pt>(8ul);
        tensor_t corr_dm{
          {2.11870860,  -0.50150667, 0.00000000, 0.00000000, 0.00000000},
          {-0.50150667, 2.11870860,  0.00000000, 0.00000000, 0.00000000},
          {0.00000000,  0.00000000,  1.33333333, 0.00000000, 0.00000000},
          {0.00000000,  0.00000000,  0.00000000, 1.33333333, 0.00000000},
          {0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.33333333}};
        atomic_basis_t corr_bs("sto-3g", 8, 0.0, 0.0, 0.0);
        corr_bs.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{1.5432896730e-01, 5.3532814230e-01,
                                           4.4463454220e-01},
                       std::vector<double>{1.3070932140e+02, 2.3808866050e+01,
                                           6.4436083130e+00});
        corr_bs.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{-9.9967229190e-02, 3.9951282610e-01,
                                           7.0011546890e-01},
                       std::vector<double>{5.0331513190e+00, 1.1695961250e+00,
                                           3.8038896000e-01});
        corr_bs.add_shell(chemist::ShellType::pure, 1,
                       std::vector<double>{1.5591627500e-01, 6.0768371860e-01,
                                           3.9195739310e-01},
                       std::vector<double>{5.0331513190e+00, 1.1695961250e+00,
                                           3.8038896000e-01});
        ao_bs_t ao_bs;
        ao_bs.add_center(corr_bs);
        aos_t corr_aos(ao_bs);

        atomic_den_t corr_rho(corr_dm, corr_aos);

        REQUIRE(rho == corr_rho);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(atom_dm_mod.run_as<atomic_den_pt>(1000ul),
                               std::out_of_range,
                               Message("Basis Set not available for Z"));
    }
}
