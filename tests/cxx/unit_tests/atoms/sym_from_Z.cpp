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
#include <simde/chemical_system/symbol_from_Z.hpp>

using sym_pt          = simde::SymbolFromZ;
using atomic_number_t = simde::type::atomic_number;

using Catch::Matchers::Message;

TEST_CASE("Symbol from Z") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& sym_mod = mm.at("Symbol from Z");
    atomic_number_t Z;

    SECTION("Hydrogen") {
        Z = 1;
        auto rv = sym_mod.run_as<sym_pt>(Z);
        REQUIRE(rv == "H");
    }

    SECTION("Oxygen") {
        Z = 8;
        auto rv = sym_mod.run_as<sym_pt>(Z);
        REQUIRE(rv == "O");
    }

    SECTION("Out of Range") {
        Z = 1000;
        REQUIRE_THROWS_MATCHES(sym_mod.run_as<sym_pt>(Z),
                               std::out_of_range,
                               Message("Symbol not available for Z"));
    }
}
