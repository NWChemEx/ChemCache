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

/*
 * This file is autogenerated by: generate_densities.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_densities.py is run.
 */

#pragma once
#include <simde/simde.hpp>

namespace chemcache {

template<typename FloatType>
simde::type::tensor sto_dash_3g_atom_density_matrix_(
  simde::type::atomic_number Z, parallelzone::runtime::RuntimeView& rt) {
    using shape_t     = tensorwrapper::shape::Smooth;
    using allocator_t = tensorwrapper::allocator::Eigen<FloatType>;
    using rank2_il_t  = typename allocator_t::rank2_il;

    allocator_t d_a(rt);
    switch(Z) {
        case(1): {
            shape_t shape{1, 1};
            auto buffer = d_a.construct(rank2_il_t{{1.00000000}});
            return simde::type::tensor(shape, std::move(buffer));
        }
        default: {
            throw std::out_of_range("Atomic Density not available for Z");
        }
    }
}

} // namespace chemcache