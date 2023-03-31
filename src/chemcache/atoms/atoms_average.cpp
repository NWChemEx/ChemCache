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
 * This file is autogenerated by: generate_atomicinfo.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_atomicinfo.py is run.
 */

#include "atoms.hpp"
#include <simde/simde.hpp>

namespace chemcache {

using atom_pt = simde::AtomFromZ;
using atom_t  = simde::type::atom;

static constexpr auto module_desc = R"(
Atoms with Abundance-Weighted Mass
---------------------------------

This module returns an atom with abundance-weighted mass.
This module was autogenerated.
)";

MODULE_CTOR(atoms_average) {
    description(module_desc);
    satisfies_property_type<atom_pt>();
}

MODULE_RUN(atoms_average) {
    const auto& [Z] = atom_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {
        case(0): {
            atom_t atom{"Ez", 0ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(1): {
            atom_t atom{"H", 1ul, 1837.4260218693814, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(2): {
            atom_t atom{"He", 2ul, 7296.297100609073, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(3): {
            atom_t atom{"Li", 3ul, 12700.97552754276, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(4): {
            atom_t atom{"Be", 4ul, 16428.204808444127, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(5): {
            atom_t atom{"B", 5ul, 19711.80464543719, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(6): {
            atom_t atom{"C", 6ul, 21893.984452257635, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(7): {
            atom_t atom{"N", 7ul, 25532.934707260847, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(8): {
            atom_t atom{"O", 8ul, 29165.122045980286, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(9): {
            atom_t atom{"F", 9ul, 34631.97038186638, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(10): {
            atom_t atom{"Ne", 10ul, 36785.3427848087, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(11): {
            atom_t atom{"Na", 11ul, 41907.785720722546, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(12): {
            atom_t atom{"Mg", 12ul, 44306.21610113965, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(13): {
            atom_t atom{"Al", 13ul, 49184.335871396164, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(14): {
            atom_t atom{"Si", 14ul, 51195.823134702325, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(15): {
            atom_t atom{"P", 15ul, 56461.71412020552, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(16): {
            atom_t atom{"S", 16ul, 58455.476530961954, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(17): {
            atom_t atom{"Cl", 17ul, 64624.13116823568, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(18): {
            atom_t atom{"Ar", 18ul, 72820.74924639802, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(19): {
            atom_t atom{"K", 19ul, 71271.84089968068, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(20): {
            atom_t atom{"Ca", 20ul, 73057.72474960299, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(21): {
            atom_t atom{"Sc", 21ul, 81949.60707950682, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(22): {
            atom_t atom{"Ti", 22ul, 87256.20316855246, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(23): {
            atom_t atom{"V", 23ul, 92860.67381934977, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(24): {
            atom_t atom{"Cr", 24ul, 94783.09201688785, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(25): {
            atom_t atom{"Mn", 25ul, 100145.9278615095, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(26): {
            atom_t atom{"Fe", 26ul, 101799.20751139225, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(27): {
            atom_t atom{"Co", 27ul, 107428.64079711946, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(28): {
            atom_t atom{"Ni", 28ul, 106991.52307546153, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(29): {
            atom_t atom{"Cu", 29ul, 115837.27174355683, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(30): {
            atom_t atom{"Zn", 30ul, 119180.44922723295, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(31): {
            atom_t atom{"Ga", 31ul, 127097.25392276482, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(32): {
            atom_t atom{"Ge", 32ul, 132396.39075212495, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(33): {
            atom_t atom{"As", 33ul, 136573.71289264012, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(34): {
            atom_t atom{"Se", 34ul, 143955.32664306843, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(35): {
            atom_t atom{"Br", 35ul, 145656.08160068557, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(36): {
            atom_t atom{"Kr", 36ul, 152754.40936591724, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(37): {
            atom_t atom{"Rb", 37ul, 155798.26856016062, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(38): {
            atom_t atom{"Sr", 38ul, 159721.48916014304, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(39): {
            atom_t atom{"Y", 39ul, 162065.43209122817, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(40): {
            atom_t atom{"Zr", 40ul, 166291.17926437902, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(41): {
            atom_t atom{"Nb", 41ul, 169357.95216689384, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(42): {
            atom_t atom{"Mo", 42ul, 174906.15025012242, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(43): {
            atom_t atom{"Tc", 43ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(44): {
            atom_t atom{"Ru", 44ul, 184239.33929942542, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(45): {
            atom_t atom{"Rh", 45ul, 187585.25111583088, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(46): {
            atom_t atom{"Pd", 46ul, 193991.79270055264, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(47): {
            atom_t atom{"Ag", 47ul, 196631.6998062559, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(48): {
            atom_t atom{"Cd", 48ul, 204918.1862867875, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(49): {
            atom_t atom{"In", 49ul, 209300.41020759306, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(50): {
            atom_t atom{"Sn", 50ul, 216395.0921958523, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(51): {
            atom_t atom{"Sb", 51ul, 221954.90207873794, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(52): {
            atom_t atom{"Te", 52ul, 232600.5708380992, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(53): {
            atom_t atom{"I", 53ul, 231332.6972092981, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(54): {
            atom_t atom{"Xe", 54ul, 239332.49801760627, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(55): {
            atom_t atom{"Cs", 55ul, 242271.81813002797, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(56): {
            atom_t atom{"Ba", 56ul, 250331.80714328878, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(57): {
            atom_t atom{"La", 57ul, 253209.1819320883, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(58): {
            atom_t atom{"Ce", 58ul, 255415.84313127832, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(59): {
            atom_t atom{"Pr", 59ul, 256858.95103025704, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(60): {
            atom_t atom{"Nd", 60ul, 262937.0810253065, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(61): {
            atom_t atom{"Pm", 61ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(62): {
            atom_t atom{"Sm", 62ul, 274089.5127838292, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(63): {
            atom_t atom{"Eu", 63ul, 277013.4259156811, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(64): {
            atom_t atom{"Gd", 64ul, 286649.214453692, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(65): {
            atom_t atom{"Tb", 65ul, 289703.1906790338, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(66): {
            atom_t atom{"Dy", 66ul, 296219.3790062, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(67): {
            atom_t atom{"Ho", 67ul, 300649.599580847, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(68): {
            atom_t atom{"Er", 68ul, 304894.5053119877, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(69): {
            atom_t atom{"Tm", 69ul, 307948.24456182634, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(70): {
            atom_t atom{"Yb", 70ul, 315441.7380930946, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(71): {
            atom_t atom{"Lu", 71ul, 318944.9651858584, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(72): {
            atom_t atom{"Hf", 72ul, 325367.3659004101, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(73): {
            atom_t atom{"Ta", 73ul, 329847.80705285165, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(74): {
            atom_t atom{"W", 74ul, 335119.8193015373, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(75): {
            atom_t atom{"Re", 75ul, 339434.5963483537, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(76): {
            atom_t atom{"Os", 76ul, 346768.07672830415, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(77): {
            atom_t atom{"Ir", 77ul, 350390.1561503677, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(78): {
            atom_t atom{"Pt", 78ul, 355616.37744028016, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(79): {
            atom_t atom{"Au", 79ul, 359048.0907948421, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(80): {
            atom_t atom{"Hg", 80ul, 365656.8472222257, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(81): {
            atom_t atom{"Tl", 81ul, 372568.3289176226, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(82): {
            atom_t atom{"Pb", 82ul, 377702.4943389824, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(83): {
            atom_t atom{"Bi", 83ul, 380947.9649997987, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(84): {
            atom_t atom{"Po", 84ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(85): {
            atom_t atom{"At", 85ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(86): {
            atom_t atom{"Rn", 86ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(87): {
            atom_t atom{"Fr", 87ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(88): {
            atom_t atom{"Ra", 88ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(89): {
            atom_t atom{"Ac", 89ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(90): {
            atom_t atom{"Th", 90ul, 422978.85169247346, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(91): {
            atom_t atom{"Pa", 91ul, 421152.64554923656, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(92): {
            atom_t atom{"U", 92ul, 433900.1594198318, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(93): {
            atom_t atom{"Np", 93ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(94): {
            atom_t atom{"Pu", 94ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(95): {
            atom_t atom{"Am", 95ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(96): {
            atom_t atom{"Cm", 96ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(97): {
            atom_t atom{"Bk", 97ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(98): {
            atom_t atom{"Cf", 98ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(99): {
            atom_t atom{"Es", 99ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(100): {
            atom_t atom{"Fm", 100ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(101): {
            atom_t atom{"Md", 101ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(102): {
            atom_t atom{"No", 102ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(103): {
            atom_t atom{"Lr", 103ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(104): {
            atom_t atom{"Rf", 104ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(105): {
            atom_t atom{"Db", 105ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(106): {
            atom_t atom{"Sg", 106ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(107): {
            atom_t atom{"Bh", 107ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(108): {
            atom_t atom{"Hs", 108ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(109): {
            atom_t atom{"Mt", 109ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(110): {
            atom_t atom{"Ds", 110ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(111): {
            atom_t atom{"Rg", 111ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(112): {
            atom_t atom{"Cn", 112ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(113): {
            atom_t atom{"Nh", 113ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(114): {
            atom_t atom{"Fl", 114ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(115): {
            atom_t atom{"Mc", 115ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(116): {
            atom_t atom{"Lv", 116ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(117): {
            atom_t atom{"Ts", 117ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        case(118): {
            atom_t atom{"Og", 118ul, 0.0, 0.0, 0.0, 0.0};
            return atom_pt::wrap_results(rv, atom);
        }
        default: {
            throw std::out_of_range("Atom not available for Z");
        }
    }
}

} // namespace chemcache
