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

#include "cc_dash_pvdz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pvdz_35() {
    // Basis Set name and origin point
    std::string name("cc-pvdz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.9840000000e-04, 1.5400000000e-03, 8.0096000000e-03,
                3.2734100000e-02, 1.0744800000e-01, 2.6889460000e-01,
                4.2044110000e-01, 2.8380410000e-01, 3.1545500000e-02,
                -7.4268000000e-03, 2.7728000000e-03, -1.3635000000e-03,
                3.8120000000e-04, -1.6150000000e-04},
      doubles_t{6.4010000000e+05, 9.5938000000e+04, 2.1833000000e+04,
                6.1819000000e+03, 2.0157000000e+03, 7.2710000000e+02,
                2.8328000000e+02, 1.1591000000e+02, 3.6124000000e+01,
                1.5532000000e+01, 4.7857000000e+00, 2.0817000000e+00,
                4.2028000000e-01, 1.6069000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.2500000000e-05, -4.8160000000e-04, -2.5466000000e-03,
                -1.0411200000e-02, -3.6517900000e-02, -9.9629500000e-02,
                -2.1413100000e-01, -1.8609110000e-01, 4.6282610000e-01,
                6.4411410000e-01, 8.2550200000e-02, -1.4969400000e-02,
                3.5288000000e-03, -1.4909000000e-03},
      doubles_t{6.4010000000e+05, 9.5938000000e+04, 2.1833000000e+04,
                6.1819000000e+03, 2.0157000000e+03, 7.2710000000e+02,
                2.8328000000e+02, 1.1591000000e+02, 3.6124000000e+01,
                1.5532000000e+01, 4.7857000000e+00, 2.0817000000e+00,
                4.2028000000e-01, 1.6069000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.4800000000e-05, 1.9190000000e-04, 1.0100000000e-03,
                4.1659000000e-03, 1.4568300000e-02, 4.0834500000e-02,
                8.9485900000e-02, 8.7278600000e-02, -2.9336440000e-01,
                -5.6671090000e-01, 5.1056580000e-01, 7.4772140000e-01,
                4.2151200000e-02, -1.0661200000e-02},
      doubles_t{6.4010000000e+05, 9.5938000000e+04, 2.1833000000e+04,
                6.1819000000e+03, 2.0157000000e+03, 7.2710000000e+02,
                2.8328000000e+02, 1.1591000000e+02, 3.6124000000e+01,
                1.5532000000e+01, 4.7857000000e+00, 2.0817000000e+00,
                4.2028000000e-01, 1.6069000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.6000000000e-06, -5.8800000000e-05, -3.0920000000e-04,
                -1.2766000000e-03, -4.4634000000e-03, -1.2557500000e-02,
                -2.7614500000e-02, -2.7394500000e-02, 9.6409400000e-02,
                1.9768710000e-01, -2.2666930000e-01, -5.2411650000e-01,
                6.8898650000e-01, 5.3443310000e-01},
      doubles_t{6.4010000000e+05, 9.5938000000e+04, 2.1833000000e+04,
                6.1819000000e+03, 2.0157000000e+03, 7.2710000000e+02,
                2.8328000000e+02, 1.1591000000e+02, 3.6124000000e+01,
                1.5532000000e+01, 4.7857000000e+00, 2.0817000000e+00,
                4.2028000000e-01, 1.6069000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6069000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.4448000000e-03, 1.2128800000e-02, 6.0807700000e-02,
                2.0093580000e-01, 4.0474190000e-01, 3.9571510000e-01,
                1.1022130000e-01, -9.0900000000e-04, 2.4832000000e-03,
                -5.7440000000e-04, 1.6910000000e-04},
      doubles_t{4.3408000000e+03, 1.0289000000e+03, 3.3202000000e+02,
                1.2516000000e+02, 5.1511000000e+01, 2.2281000000e+01,
                9.3417000000e+00, 4.0132000000e+00, 1.7002000000e+00,
                4.7194000000e-01, 1.4421000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-5.8190000000e-04, -4.9065000000e-03, -2.5251400000e-02,
                -8.6944500000e-02, -1.8934220000e-01, -1.7108820000e-01,
                2.3687550000e-01, 5.8984000000e-01, 3.1719440000e-01,
                1.7983300000e-02, -1.4683000000e-03},
      doubles_t{4.3408000000e+03, 1.0289000000e+03, 3.3202000000e+02,
                1.2516000000e+02, 5.1511000000e+01, 2.2281000000e+01,
                9.3417000000e+00, 4.0132000000e+00, 1.7002000000e+00,
                4.7194000000e-01, 1.4421000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.5180000000e-04, 1.2563000000e-03, 6.6224000000e-03,
                2.2381600000e-02, 5.0971700000e-02, 4.1400900000e-02,
                -7.0397000000e-02, -2.2325400000e-01, -5.6417900000e-02,
                5.8080790000e-01, 5.5081320000e-01},
      doubles_t{4.3408000000e+03, 1.0289000000e+03, 3.3202000000e+02,
                1.2516000000e+02, 5.1511000000e+01, 2.2281000000e+01,
                9.3417000000e+00, 4.0132000000e+00, 1.7002000000e+00,
                4.7194000000e-01, 1.4421000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.4421000000e-01}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{2.2658300000e-02, 1.3458950000e-01, 3.6471810000e-01,
                           4.9041960000e-01, 2.7138850000e-01},
                 doubles_t{1.0483000000e+02, 3.0272000000e+01, 1.0649000000e+01,
                           3.8696000000e+00, 1.3239000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0980000000e-01}));
    return abs_t(name, 35, r0, shells.begin(), shells.end());
} // cc_dash_pvdz_35

} // namespace chemcache
