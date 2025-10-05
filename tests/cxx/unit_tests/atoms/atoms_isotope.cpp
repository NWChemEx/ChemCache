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

#include "../catch.hpp"
#include "chemcache/chemcache.hpp"
#include <simde/chemical_system/atom.hpp>
#include <simde/types.hpp>

using z_t        = simde::type::atomic_number;
using input_t    = std::pair<z_t, z_t>;
using isotope_pt = simde::Atom<input_t>;
using atom_t     = simde::type::atom;

using Catch::Matchers::Message;

TEST_CASE("Atom Isotope") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& atom_mod = mm.at("Atom Isotope");

    SECTION("Hydrogen-2") {
        input_t in{1, 2};
        auto rv = atom_mod.run_as<isotope_pt>(in);
        atom_t corr{"H", 1, 3671.4829413173247, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Oxygen-18") {
        input_t in{8, 18};
        auto rv = atom_mod.run_as<isotope_pt>(in);
        atom_t corr{"O", 8, 32810.46081966976, 0.0, 0.0, 0.0};
        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        SECTION("Z Out of Range") {
            input_t in{1000, 1};
            REQUIRE_THROWS_MATCHES(atom_mod.run_as<isotope_pt>(in),
                                   std::out_of_range,
                                   Message("Isotopes not available for Z: " +
                                           std::to_string(in.first)));
        }
        SECTION("N Out of Range") {
            input_t in{1, 10000};
            auto z_string = std::to_string(in.first);
            auto n_string = std::to_string(in.second);
            REQUIRE_THROWS_MATCHES(
              atom_mod.run_as<isotope_pt>(in), std::out_of_range,
              Message("Isotope not available for Z and mass number: (" +
                      z_string + ", " + n_string + ")"));
        }
    }
}
