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

using Catch::Matchers::Message;

TEST_CASE("Atomic Basis Set") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_bs_mod = mm.at("sto-3g atomic basis");

    SECTION("sto-3g Hydrogen") {
        auto [rv] = atom_bs_mod.run_as<atomic_basis_pt>(1ul);
        atomic_basis_t corr("sto-3g", 1, 0.0, 0.0, 0.0);
        corr.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{1.5432896730e-01, 5.3532814230e-01,
                                           4.4463454220e-01},
                       std::vector<double>{3.4252509140e+00, 6.2391372980e-01,
                                           1.6885540400e-01});
        REQUIRE(rv == corr);
    }

    SECTION("sto-3g Oxygen") {
        auto [rv] = atom_bs_mod.run_as<atomic_basis_pt>(8ul);
        atomic_basis_t corr("sto-3g", 8, 0.0, 0.0, 0.0);
        corr.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{1.5432896730e-01, 5.3532814230e-01,
                                           4.4463454220e-01},
                       std::vector<double>{1.3070932140e+02, 2.3808866050e+01,
                                           6.4436083130e+00});
        corr.add_shell(chemist::ShellType::pure, 0,
                       std::vector<double>{-9.9967229190e-02, 3.9951282610e-01,
                                           7.0011546890e-01},
                       std::vector<double>{5.0331513190e+00, 1.1695961250e+00,
                                           3.8038896000e-01});
        corr.add_shell(chemist::ShellType::pure, 1,
                       std::vector<double>{1.5591627500e-01, 6.0768371860e-01,
                                           3.9195739310e-01},
                       std::vector<double>{5.0331513190e+00, 1.1695961250e+00,
                                           3.8038896000e-01});
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(atom_bs_mod.run_as<atomic_basis_pt>(1000ul),
                               std::out_of_range,
                               Message("Basis Set not available for Z"));
    }
}