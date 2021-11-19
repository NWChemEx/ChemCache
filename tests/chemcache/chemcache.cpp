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