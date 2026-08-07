/*
 * Copyright 2025 NWChemEx-Project
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

#include "density.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using sad_rho_pt      = simde::InitialDensity;
using atomic_den_pt   = simde::AtomicDensityMatrixFromZ;
using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using aobs_t          = simde::type::ao_basis_set;
using cmos_t          = simde::type::cmos;
using aos_t           = simde::type::aos;
using density_t       = simde::type::decomposable_e_density;
using tensor_t        = simde::type::tensor;

namespace {

// Pull out nuclear-nuclear interaction term, if there is one.
struct GrabNuclear : chemist::qm_operator::OperatorVisitor {
    using V_en_type = simde::type::V_en_type;

    GrabNuclear() : chemist::qm_operator::OperatorVisitor(false) {}

    void run(const V_en_type& V_en) { m_pv = &V_en; }

    const V_en_type* m_pv = nullptr;
};

} // namespace

MODULE_CTOR(sad_density) {
    satisfies_property_type<sad_rho_pt>();
    add_submodule<atomic_den_pt>("Atomic Density");
    add_submodule<atomic_basis_pt>("Atomic Basis");
}

MODULE_RUN(sad_density) {
    const auto& [H]   = sad_rho_pt::unwrap_inputs(inputs);
    auto& atom_dm_mod = submods.at("Atomic Density");
    auto& atom_bs_mod = submods.at("Atomic Basis");

    // get atomic info from H
    GrabNuclear visitor;
    H.visit(visitor);
    if(visitor.m_pv == nullptr)
        throw std::runtime_error("No nuclear terms for SAD density generation");
    const auto nuclei = visitor.m_pv->get_rhs_particle();

    // Aggregate basis sets and tensors from centers
    aobs_t aobs;
    std::vector<simde::type::tensor> tensors(nuclei.size());
    for(std::size_t ni = 0; ni < nuclei.size(); ++ni) {
        auto Z      = nuclei[ni].Z();
        tensors[ni] = atom_dm_mod.run_as<atomic_den_pt>(Z);
        auto ci     = atom_bs_mod.run_as<atomic_basis_pt>(Z);
        for(auto i : {0, 1, 2}) ci.center().coord(i) = nuclei[ni].coord(i);
        aobs.add_center(ci);
    }

    // Stack density matrices and scale
    auto rho    = tensorwrapper::utilities::block_diagonal_matrix(tensors);
    rho("i, j") = rho("i, j") * 0.5;

    aos_t aos(aobs);
    cmos_t cmos(tensor_t{}, aos, tensor_t{});
    density_t d(rho, cmos);

    auto rv = results();
    return sad_rho_pt::wrap_results(rv, d);
}

} // namespace chemcache
