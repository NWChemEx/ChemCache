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

#include <catch2/catch.hpp>
#include <chemcache/chemcache.hpp>
#include <chemist/chemist.hpp>

using namespace chemcache;

TEST_CASE("Load into managers") {
    // These tests are just going to make sure that the code compiles and
    // uses the managers correctly
    SECTION("load_basis_sets") {
        chemist::BasisSetManager bsm;

        REQUIRE_NOTHROW(load_basis_sets(bsm));
    }

    SECTION("load_molecules") {
        chemist::MoleculeManager mm;
        chemist::PeriodicTable pt;

        load_elements(pt);
        
        REQUIRE_NOTHROW(load_molecules(mm, pt));
    }

    SECTION("load_elements") {
        chemist::PeriodicTable pt;

        REQUIRE_NOTHROW(load_elements(pt));
    }
}

TEST_CASE("Create default managers") {
    // These tests are just going to make sure that the code compiles and
    // uses the managers correctly
    SECTION("nwx_basis_set_manager") {
        REQUIRE_NOTHROW(nwx_basis_set_manager());
    }

    SECTION("nwx_molecule_manager") {
        REQUIRE_NOTHROW(nwx_molecule_manager());
    }

    SECTION("nwx_periodic_table") {
        REQUIRE_NOTHROW(nwx_periodic_table());

    }
}