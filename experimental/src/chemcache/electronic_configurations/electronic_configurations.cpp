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
 * This file is autogenerated by: generate_elec_configs.py 
 * 
 * NOTE: Any modifications made in this file will be lost next time 
 *       generate_elec_configs.py is run.
 */

#include "electronic_configurations.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using elec_config_pt = simde::ElecConfigFromZ;
using elec_config_t  = std::vector<simde::type::size>;

static constexpr auto module_desc = R"(
Electronic Configurations
---------------------------------

This module returns the electron configuration associated with the atomic 
number. This module was autogenerated.
)";

MODULE_CTOR(elec_configs) {
    description(module_desc);
    satisfies_property_type<elec_config_pt>();
}

MODULE_RUN(elec_configs) {
    const auto& [Z] = elec_config_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {
        case(0): {
            elec_config_t elec_config{0, 0, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(1): {
            elec_config_t elec_config{1, 0, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(2): {
            elec_config_t elec_config{2, 0, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(3): {
            elec_config_t elec_config{3, 0, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(4): {
            elec_config_t elec_config{4, 0, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(5): {
            elec_config_t elec_config{4, 1, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(6): {
            elec_config_t elec_config{4, 2, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(7): {
            elec_config_t elec_config{4, 3, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(8): {
            elec_config_t elec_config{4, 4, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(9): {
            elec_config_t elec_config{4, 5, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(10): {
            elec_config_t elec_config{4, 6, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(11): {
            elec_config_t elec_config{5, 6, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(12): {
            elec_config_t elec_config{6, 6, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(13): {
            elec_config_t elec_config{6, 7, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(14): {
            elec_config_t elec_config{6, 8, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(15): {
            elec_config_t elec_config{6, 9, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(16): {
            elec_config_t elec_config{6, 10, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(17): {
            elec_config_t elec_config{6, 11, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(18): {
            elec_config_t elec_config{6, 12, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(19): {
            elec_config_t elec_config{7, 12, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(20): {
            elec_config_t elec_config{8, 12, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(21): {
            elec_config_t elec_config{8, 13, 0, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(22): {
            elec_config_t elec_config{8, 12, 2, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(23): {
            elec_config_t elec_config{8, 12, 3, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(24): {
            elec_config_t elec_config{8, 12, 4, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(25): {
            elec_config_t elec_config{6, 12, 7, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(26): {
            elec_config_t elec_config{6, 12, 8, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(27): {
            elec_config_t elec_config{6, 12, 9, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(28): {
            elec_config_t elec_config{6, 12, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(29): {
            elec_config_t elec_config{7, 12, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(30): {
            elec_config_t elec_config{8, 12, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(31): {
            elec_config_t elec_config{8, 13, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(32): {
            elec_config_t elec_config{8, 14, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(33): {
            elec_config_t elec_config{8, 15, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(34): {
            elec_config_t elec_config{8, 16, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(35): {
            elec_config_t elec_config{8, 17, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(36): {
            elec_config_t elec_config{8, 18, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(37): {
            elec_config_t elec_config{9, 18, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(38): {
            elec_config_t elec_config{10, 18, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(39): {
            elec_config_t elec_config{10, 19, 10, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(40): {
            elec_config_t elec_config{10, 18, 12, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(41): {
            elec_config_t elec_config{10, 18, 13, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(42): {
            elec_config_t elec_config{8, 18, 16, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(43): {
            elec_config_t elec_config{8, 18, 17, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(44): {
            elec_config_t elec_config{8, 18, 18, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(45): {
            elec_config_t elec_config{8, 18, 19, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(46): {
            elec_config_t elec_config{8, 18, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(47): {
            elec_config_t elec_config{9, 18, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(48): {
            elec_config_t elec_config{10, 18, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(49): {
            elec_config_t elec_config{10, 19, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(50): {
            elec_config_t elec_config{10, 20, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(51): {
            elec_config_t elec_config{10, 21, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(52): {
            elec_config_t elec_config{10, 22, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(53): {
            elec_config_t elec_config{10, 23, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(54): {
            elec_config_t elec_config{10, 24, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(55): {
            elec_config_t elec_config{11, 24, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(56): {
            elec_config_t elec_config{12, 24, 20, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(57): {
            elec_config_t elec_config{12, 24, 21, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(58): {
            elec_config_t elec_config{12, 24, 22, 0};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(59): {
            elec_config_t elec_config{12, 24, 21, 2};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(60): {
            elec_config_t elec_config{12, 24, 20, 4};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(61): {
            elec_config_t elec_config{12, 24, 20, 5};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(62): {
            elec_config_t elec_config{12, 24, 20, 6};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(63): {
            elec_config_t elec_config{12, 24, 20, 7};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(64): {
            elec_config_t elec_config{11, 24, 20, 9};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(65): {
            elec_config_t elec_config{10, 24, 20, 11};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(66): {
            elec_config_t elec_config{10, 24, 20, 12};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(67): {
            elec_config_t elec_config{10, 24, 20, 13};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(68): {
            elec_config_t elec_config{10, 24, 20, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(69): {
            elec_config_t elec_config{11, 24, 20, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(70): {
            elec_config_t elec_config{12, 24, 20, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(71): {
            elec_config_t elec_config{12, 25, 20, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(72): {
            elec_config_t elec_config{12, 24, 22, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(73): {
            elec_config_t elec_config{12, 24, 23, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(74): {
            elec_config_t elec_config{10, 24, 26, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(75): {
            elec_config_t elec_config{10, 24, 27, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(76): {
            elec_config_t elec_config{10, 24, 28, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(77): {
            elec_config_t elec_config{10, 24, 29, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(78): {
            elec_config_t elec_config{10, 24, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(79): {
            elec_config_t elec_config{11, 24, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(80): {
            elec_config_t elec_config{12, 24, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(81): {
            elec_config_t elec_config{12, 25, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(82): {
            elec_config_t elec_config{12, 26, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(83): {
            elec_config_t elec_config{12, 27, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(84): {
            elec_config_t elec_config{12, 28, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(85): {
            elec_config_t elec_config{12, 29, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(86): {
            elec_config_t elec_config{12, 30, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(87): {
            elec_config_t elec_config{13, 30, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(88): {
            elec_config_t elec_config{14, 30, 30, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(89): {
            elec_config_t elec_config{14, 30, 31, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(90): {
            elec_config_t elec_config{14, 30, 32, 14};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(91): {
            elec_config_t elec_config{14, 30, 30, 17};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(92): {
            elec_config_t elec_config{14, 30, 30, 18};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(93): {
            elec_config_t elec_config{14, 30, 30, 19};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(94): {
            elec_config_t elec_config{13, 30, 30, 21};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(95): {
            elec_config_t elec_config{12, 30, 30, 23};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(96): {
            elec_config_t elec_config{12, 30, 30, 24};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(97): {
            elec_config_t elec_config{12, 30, 30, 25};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(98): {
            elec_config_t elec_config{12, 30, 30, 26};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(99): {
            elec_config_t elec_config{12, 30, 30, 27};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(100): {
            elec_config_t elec_config{12, 30, 30, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(101): {
            elec_config_t elec_config{13, 30, 30, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(102): {
            elec_config_t elec_config{14, 30, 30, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(103): {
            elec_config_t elec_config{14, 30, 31, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(104): {
            elec_config_t elec_config{14, 30, 32, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(105): {
            elec_config_t elec_config{14, 30, 33, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(106): {
            elec_config_t elec_config{12, 30, 36, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(107): {
            elec_config_t elec_config{12, 30, 37, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(108): {
            elec_config_t elec_config{12, 30, 38, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(109): {
            elec_config_t elec_config{12, 30, 39, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(110): {
            elec_config_t elec_config{12, 30, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(111): {
            elec_config_t elec_config{13, 30, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(112): {
            elec_config_t elec_config{14, 30, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(113): {
            elec_config_t elec_config{14, 31, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(114): {
            elec_config_t elec_config{14, 32, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(115): {
            elec_config_t elec_config{14, 33, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(116): {
            elec_config_t elec_config{14, 34, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(117): {
            elec_config_t elec_config{14, 35, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        case(118): {
            elec_config_t elec_config{14, 36, 40, 28};
            return elec_config_pt::wrap_results(rv, elec_config);
        }
        default: {
            throw std::out_of_range("Atomic Density not available for Z");
        }
    }
}

} // namespace chemcache