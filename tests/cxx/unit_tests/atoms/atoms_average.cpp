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

using atom_pt = simde::AtomFromZ;
using atom_t  = simde::type::atom;

using Catch::Matchers::Message;

TEST_CASE("Atom") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_mod = mm.at("Atom");

    SECTION("Hydrogen") {
        auto rv = atom_mod.run_as<atom_pt>(1ul);
        atom_t corr{"H", 1ul, 1837.4260218693814, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Oxygen") {
        auto rv = atom_mod.run_as<atom_pt>(8ul);
        atom_t corr{"O", 8ul, 29165.122045980286, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(atom_mod.run_as<atom_pt>(1000ul),
                               std::out_of_range,
                               Message("Atom not available for Z"));
    }
}
