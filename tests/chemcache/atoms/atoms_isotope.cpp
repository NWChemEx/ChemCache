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

using size_t     = simde::type::size;
using input_t    = std::pair<size_t, size_t>;
using isotope_pt = simde::Atom<input_t>;
using atom_t     = simde::type::atom;

using Catch::Matchers::Message;

TEST_CASE("Atom Isotope") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_mod = mm.at("Atom Isotope");

    SECTION("Hydrogen-2") {
        input_t in{1ul, 2ul};
        auto rv = atom_mod.run_as<isotope_pt>(in);
        atom_t corr{"H", 1ul, 3671.4829413173247, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Oxygen-18") {
        input_t in{8ul, 18ul};
        auto rv = atom_mod.run_as<isotope_pt>(in);
        atom_t corr{"O", 8ul, 32810.46081966976, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        SECTION("Z Out of Range") {
            input_t in{1000ul, 1ul};
            REQUIRE_THROWS_MATCHES(
              atom_mod.run_as<isotope_pt>(in), std::out_of_range,
              Message("Isotopes not available for Z"));
        }
        SECTION("N Out of Range") {
            input_t in{1ul, 10000ul};
            REQUIRE_THROWS_MATCHES(
              atom_mod.run_as<isotope_pt>(in), std::out_of_range,
              Message("Isotope not available for Z and mass number"));
        }
    }
}
