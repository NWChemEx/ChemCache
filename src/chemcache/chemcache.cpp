/**
 * @file chemcache.cpp
 *
 * @brief Implementation of chemcache.hpp API.
 */

#include "chemcache/chemcache.hpp"

namespace chemcache {

libchemist::BasisSetManager nwx_basis_set_manager() {
    // TODO: Update with chemcache defaults
    return libchemist::BasisSetManager();
}

libchemist::MoleculeManager nwx_molecule_manager() {
    // TODO: Update with chemcache defaults
    return libchemist::MoleculeManager();
}

libchemist::PeriodicTable nwx_periodic_table() {
    // TODO: Update with chemcache defaults
    return libchemist::PeriodicTable();
}

} // namespace chemcache