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

using sym_pt = simde::SymbolFromZ;

using Catch::Matchers::Message;

TEST_CASE("Symbol from Z") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& sym_mod = mm.at("Symbol from Z");

    SECTION("Hydrogen") {
        auto rv = sym_mod.run_as<sym_pt>(1ul);
        REQUIRE(rv == "H");
    }

    SECTION("Oxygen") {
        auto rv = sym_mod.run_as<sym_pt>(8ul);
        REQUIRE(rv == "O");
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(sym_mod.run_as<sym_pt>(1000ul),
                               std::out_of_range,
                               Message("Symbol not available for Z"));
    }
}
