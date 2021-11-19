/**
 * @file chemcache.cpp
 *
 * @brief Implementation of chemcache.hpp API
 */

#include "chemcache/chemcache.hpp"

namespace chemcache {

chemist::BasisSetManager nwx_basis_set_manager() {
    chemist::BasisSetManager bsm;
    load_basis_sets(bsm);

    return bsm;
}

chemist::MoleculeManager nwx_molecule_manager() {
    chemist::MoleculeManager mm;
    chemist::PeriodicTable pt;

    load_elements(pt);
    load_molecules(mm, pt);

    return mm;
}

chemist::PeriodicTable nwx_periodic_table() {
    chemist::PeriodicTable pt;
    load_elements(pt);

    return pt;
}

} // namespace chemcache