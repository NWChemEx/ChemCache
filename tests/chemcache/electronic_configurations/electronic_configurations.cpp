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

using elec_config_pt = simde::ElecConfigFromZ;
using elec_config_t  = std::vector<simde::type::size>;

using Catch::Matchers::Message;

TEST_CASE("Electronic Configurations") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& elec_config_mod = mm.at("Elec Config from Z");

    SECTION("Hydrogen") {
        auto [rv] = elec_config_mod.run_as<elec_config_pt>(1ul);
        elec_config_t corr{1, 0, 0, 0};
        REQUIRE(rv == corr);
    }

    SECTION("Oxygen") {
        auto [rv] = elec_config_mod.run_as<elec_config_pt>(8ul);
        elec_config_t corr{4, 4, 0, 0};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        REQUIRE_THROWS_MATCHES(elec_config_mod.run_as<elec_config_pt>(1000ul),
                               std::out_of_range,
                               Message("Atomic Density not available for Z"));
    }
}