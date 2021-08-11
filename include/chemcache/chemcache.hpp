/**
 * @file chemcache.hpp
 *
 * This is the main header of the chemcache library, defining the public API
 * of the library. This file should NOT be included in any other chemcache
 * header files, source files, or tests (except tests/chemcache.cpp, which
 * tests the functions defined in this header file).
 */

#pragma once

#include <libchemist/libchemist.hpp>

/**
 * @brief The primary namespace for the chemcache library.
 *
 * This namespace contains all functionality of the chemcache library available
 * to an end-user of the library.
 */
namespace chemcache {

libchemist::BasisSetManager nwx_basis_set_manager();
libchemist::MoleculeManager nwx_molecule_manager();
libchemist::PeriodicTable nwx_periodic_table();

} // namespace chemcache