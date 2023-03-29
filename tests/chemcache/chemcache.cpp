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

    SECTION("load_elec_configs") {
        chemist::PeriodicTable pt;

        load_elements(pt);

        REQUIRE_NOTHROW(load_elec_configs(pt));
    }

    SECTION("load_atom_dm, non-default case") {
        chemist::PeriodicTable pt;

        pt.insert(8,
                  chemist::Atom("O", 8ul, 29165.122045980286, 0.0, 0.0, 0.0));

        REQUIRE_NOTHROW(load_atom_dm(6, "sto-3g", pt));
        REQUIRE_NOTHROW(load_atom_dm(7, "3-21g", pt));
        REQUIRE_NOTHROW(load_atom_dm(8, "6-31g", pt));

        chemist::PeriodicTable::atom_dm_t corr1 = {
          2.13147782,  -0.52937893, 0.00000000, 0.00000000, 0.00000000,
          -0.52937893, 2.13147782,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  0.66666667, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 0.66666667, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 0.66666667};
        chemist::PeriodicTable::atom_dm_t corr2 = {
          2.06347599,  0.04789290, 0.00000000, 0.00000000, 0.00000000,
          -0.45356235, 0.00000000, 0.00000000, 0.00000000, 0.04789290,
          0.14535500,  0.00000000, 0.00000000, 0.00000000, 0.40960696,
          0.00000000,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.30390186,  0.00000000, 0.00000000, 0.00000000, 0.32912464,
          0.00000000,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.30390186,  0.00000000, 0.00000000, 0.00000000, 0.32912464,
          0.00000000,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.30390186,  0.00000000, 0.00000000, 0.00000000, 0.32912464,
          -0.45356235, 0.40960696, 0.00000000, 0.00000000, 0.00000000,
          1.32340885,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.32912464, 0.00000000, 0.00000000, 0.00000000,
          0.35644083,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.32912464, 0.00000000, 0.00000000, 0.00000000,
          0.35644083,  0.00000000, 0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.32912464, 0.00000000, 0.00000000, 0.00000000,
          0.35644083};
        chemist::PeriodicTable::atom_dm_t corr3 = {
          2.09261376,  -0.21864494, 0.00000000, 0.00000000, 0.00000000,
          -0.25894505, 0.00000000,  0.00000000, 0.00000000, -0.21864494,
          0.61469296,  0.00000000,  0.00000000, 0.00000000, 0.59376175,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.64037418,  0.00000000,  0.00000000, 0.00000000, 0.41836429,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.64037418,  0.00000000,  0.00000000, 0.00000000, 0.41836429,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.64037418,  0.00000000,  0.00000000, 0.00000000, 0.41836429,
          -0.25894505, 0.59376175,  0.00000000, 0.00000000, 0.00000000,
          0.57467468,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.41836429,  0.00000000, 0.00000000, 0.00000000,
          0.27332251,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.41836429,  0.00000000, 0.00000000, 0.00000000,
          0.27332251,  0.00000000,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.41836429,  0.00000000, 0.00000000, 0.00000000,
          0.27332251};

        REQUIRE(corr1 == pt.get_atom_dm(6, "sto-3g"));
        REQUIRE(corr2 == pt.get_atom_dm(7, "3-21g"));
        REQUIRE(
          corr3 ==
          pt.get_atom_dm("O", "6-31g")); // failed if O is not explicitly loaded
    }

    SECTION("load_atom_dm, default case") {
        chemist::PeriodicTable pt;

        REQUIRE_NOTHROW(load_atom_dm(pt));

        chemist::PeriodicTable::atom_dm_t corr1 = {
          2.08820832,  -0.42918219, 0.00000000, 0.00000000, 0.00000000,
          -0.42918219, 2.08820832,  0.00000000, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  1.66666667, 0.00000000, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 1.66666667, 0.00000000,
          0.00000000,  0.00000000,  0.00000000, 0.00000000, 1.66666667};

        REQUIRE(corr1 == pt.get_atom_dm(9, "mini"));
    }
}

TEST_CASE("Create default managers") {
    // These tests are just going to make sure that the code compiles and
    // uses the managers correctly
    SECTION("nwx_basis_set_manager") {
        REQUIRE_NOTHROW(nwx_basis_set_manager());
    }

    SECTION("nwx_molecule_manager") { REQUIRE_NOTHROW(nwx_molecule_manager()); }

    SECTION("nwx_periodic_table") { REQUIRE_NOTHROW(nwx_periodic_table()); }
}
