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

using molecule_pt = simde::MoleculeFromString;
using atom_pt     = simde::AtomFromZ;
using molecule_t  = simde::type::molecule;
using atom_t      = simde::type::atom;

using Catch::Matchers::Message;

TEST_CASE("NWX Molecules") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    auto& mol_mod = mm.at("NWX Molecules");

    SECTION("Water") {
        double o_mass = 29165.122045980286;
        double h_mass = 1837.4260218693814;

        auto atom_mod = pluginplay::make_lambda<atom_pt>([=](auto&& Z) {
            REQUIRE(((Z == 8) || (Z == 1)));
            if(Z == 8) return atom_t{"O", 8ul, o_mass, 0.0, 0.0, 0.0};
            return atom_t{"H", 1ul, h_mass, 0.0, 0.0, 0.0};
        });
        mol_mod.change_submod("Atoms", atom_mod);

        std::string name{"water"};
        auto rv = mol_mod.run_as<molecule_pt>(name);

        double o_y = -0.1432223429807816;
        double h_x = 1.6380335020342418;
        double h_y = 1.1365568803584036;
        atom_t o{"O", 8ul, o_mass, 0.0, o_y, 0.0};
        atom_t h1{"H", 1ul, h_mass, h_x, h_y, 0.0};
        atom_t h2{"H", 1ul, h_mass, -h_x, h_y, 0.0};
        molecule_t corr{o, h1, h2};

        REQUIRE(rv == corr);
    }

    SECTION("Out of Range") {
        std::string name{"Nothing"};
        REQUIRE_THROWS_MATCHES(mol_mod.run_as<molecule_pt>(name),
                               std::out_of_range,
                               Message("No molecule found for name: Nothing"));
    }
}
