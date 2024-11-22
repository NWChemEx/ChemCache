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

/*
 * This file is autogenerated by: generate_basis.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis.py is run.
 */

#include "../bases.hpp"
#include "cc_dash_pv5z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
cc-pv5z atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR(cc_dash_pv5z_atom_basis) {
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}

MODULE_RUN(cc_dash_pv5z_atom_basis) {
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {
        case(1): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_1());
        }
        case(2): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_2());
        }
        case(3): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_3());
        }
        case(4): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_4());
        }
        case(5): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_5());
        }
        case(6): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_6());
        }
        case(7): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_7());
        }
        case(8): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_8());
        }
        case(9): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_9());
        }
        case(10): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_10());
        }
        case(11): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_11());
        }
        case(12): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_12());
        }
        case(13): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_13());
        }
        case(14): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_14());
        }
        case(15): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_15());
        }
        case(16): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_16());
        }
        case(17): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_17());
        }
        case(18): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_18());
        }
        case(20): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_20());
        }
        case(21): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_21());
        }
        case(22): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_22());
        }
        case(23): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_23());
        }
        case(24): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_24());
        }
        case(25): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_25());
        }
        case(26): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_26());
        }
        case(27): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_27());
        }
        case(28): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_28());
        }
        case(29): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_29());
        }
        case(30): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_30());
        }
        case(31): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_31());
        }
        case(32): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_32());
        }
        case(33): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_33());
        }
        case(34): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_34());
        }
        case(35): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_35());
        }
        case(36): {
            return atomic_basis_pt::wrap_results(rv, cc_dash_pv5z_36());
        }
        default: {
            throw std::out_of_range("Basis Set not available for Z");
        }
    }
}

} // namespace chemcache
