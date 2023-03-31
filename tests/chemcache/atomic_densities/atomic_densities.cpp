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

using atomic_den_pt = simde::AtomDenFromZ;
using atomic_den_t  = std::vector<double>;

using Catch::Matchers::Message;

TEST_CASE("Atomic Densities") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_dm_mod = mm.at("sto-3g atomic dm");

    SECTION("sto-3g Hydrogen") {
        auto [rv] = atom_dm_mod.run_as<atomic_den_pt>(1ul);
        atomic_den_t corr{1.00000000};
        REQUIRE(rv == corr);
    }

    SECTION("sto-3g Oxygen") {
        auto [rv] = atom_dm_mod.run_as<atomic_den_pt>(8ul);
        atomic_den_t corr{
          2.11870860,  -0.50150667, 0.00000000, 0.00000000, 0.00000000,
          -0.50150667, 2.11870860,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  1.33333333, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 1.33333333, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.33333333};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(atom_dm_mod.run_as<atomic_den_pt>(1000ul),
                               std::out_of_range,
                               Message("Atomic Density not available for Z"));
    }
}