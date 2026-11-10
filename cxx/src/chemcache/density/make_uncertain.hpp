/*
 * Copyright 2026 NWChemEx-Project
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

#pragma once
#include <simde/simde.hpp>
#include <vector>

namespace chemcache {

/** @brief Converts a tensor of doubles to the kind of values @p factory makes.
 *
 *  The values of @p t are exact, so each becomes an uncertain value with a
 *  radius of zero. This is needed because WTF does not convert doubles to
 *  uncertain values in mixed arithmetic.
 *
 *  @param[in] t The tensor to convert. Must hold doubles.
 *  @param[in] factory The factory making the values of the result.
 *
 *  @return A tensor with the values of @p t, of the kind @p factory makes.
 *
 *  @throw std::bad_alloc if there is a problem allocating the result. Strong
 *                        throw guarantee.
 */
inline simde::type::tensor make_uncertain(const simde::type::tensor& t,
                                          const simde::UQFactory& factory) {
    if(factory.kind() == simde::UQKind::none) return t;

    using tensorwrapper::buffer::get_raw_data;
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(t.buffer());
    const auto values  = get_raw_data<double>(buffer);
    auto shape         = t.buffer().layout().shape().as_smooth().make_smooth();

    return simde::visit_uq_kind(factory.kind(), [&](auto type) {
        using float_type = typename decltype(type)::type;
        std::vector<float_type> data;
        data.reserve(values.size());
        for(const auto x : values) {
            auto value = factory(x, 0.0);
            data.push_back(wtf::fp::float_cast<float_type>(value));
        }
        tensorwrapper::buffer::Contiguous uq_buffer(std::move(data), shape);
        return simde::type::tensor(shape, std::move(uq_buffer));
    });
}

} // namespace chemcache
