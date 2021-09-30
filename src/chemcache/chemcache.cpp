/**
 * @file chemcache.cpp
 *
 * @brief Implementation of chemcache.hpp API
 */

#include "chemcache/chemcache.hpp"

namespace chemcache {

libchemist::BasisSetManager nwx_basis_set_manager() {
    libchemist::BasisSetManager bsm;
    load_basis_sets(bsm);

    return bsm;
}

libchemist::MoleculeManager nwx_molecule_manager() {
    libchemist::MoleculeManager mm;
    libchemist::PeriodicTable pt;

    load_elements(pt);
    load_molecules(mm, pt);

    return mm;
}

libchemist::PeriodicTable nwx_periodic_table() {
    libchemist::PeriodicTable pt;
    load_elements(pt);

    return pt;
}

} // namespace chemcache