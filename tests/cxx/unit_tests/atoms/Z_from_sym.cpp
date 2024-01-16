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
#include <simde/atoms/Z_from_symbol.hpp>

using z_pt = simde::ZFromSymbol;

using Catch::Matchers::Message;

TEST_CASE("Z from Symbol") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& z_mod = mm.at("Z From Symbol");

    SECTION("sto-3g Hydrogen") {
        auto rv = z_mod.run_as<z_pt>("H");
        REQUIRE(rv == 1ul);
    }

    SECTION("sto-3g Oxygen") {
        auto rv = z_mod.run_as<z_pt>("O");
        REQUIRE(rv == 8ul);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(z_mod.run_as<z_pt>("Not a symbol"),
                               std::out_of_range,
                               Message("Z not available for Symbol"));
    }
}
