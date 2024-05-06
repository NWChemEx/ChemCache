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
#include <simde/chemical_system/atom.hpp>
#include <simde/types.hpp>

using atom_pt         = simde::AtomFromZ;
using atom_t          = simde::type::atom;
using atomic_number_t = simde::type::atomic_number;

using Catch::Matchers::Message;

TEST_CASE("Atom") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_mod = mm.at("Atom");
    atomic_number_t Z;

    SECTION("Hydrogen") {
        Z = 1;
        auto rv = atom_mod.run_as<atom_pt>(Z);
        atom_t corr{"H", 1, 1837.4260218693814, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Oxygen") {
        Z = 8;
        auto rv = atom_mod.run_as<atom_pt>(Z);
        atom_t corr{"O", 8, 29165.122045980286, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        Z = 1000;
        REQUIRE_THROWS_MATCHES(atom_mod.run_as<atom_pt>(Z),
                               std::out_of_range,
                               Message("Atom not available for Z"));
    }
}
