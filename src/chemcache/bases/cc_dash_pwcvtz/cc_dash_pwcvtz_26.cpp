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

#include "cc_dash_pwcvtz.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t cc_dash_pwcvtz_26() {
    // Basis Set name and origin point
    std::string name("cc-pwcvtz");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{8.0488030000e-06,  6.2583060000e-05, 3.2902390000e-04,
                1.3873550000e-03,  5.0232560000e-03, 1.6101400000e-02,
                4.5900340000e-02,  1.1361540000e-01, 2.2838690000e-01,
                3.2211590000e-01,  2.3836610000e-01, 7.4046670000e-02,
                9.2141970000e-02,  9.3397900000e-02, 1.5739650000e-02,
                -4.1866820000e-04, 5.3763180000e-05, -3.8166540000e-05,
                4.3196030000e-05,  -3.4010190000e-06},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.1559540000e-06, -3.2314010000e-05, -1.6995250000e-04,
                -7.1713690000e-04, -2.6036250000e-03, -8.3991090000e-03,
                -2.4341090000e-02, -6.2519480000e-02, -1.3659290000e-01,
                -2.3127070000e-01, -2.3837340000e-01, 3.1238370000e-02,
                5.0868180000e-01,  4.9876950000e-01,  9.0335520000e-02,
                -6.0053370000e-03, 2.3124540000e-04,  -5.6436800000e-04,
                4.9922600000e-04,  -1.0152930000e-04},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.5321780000e-07,  7.4146050000e-06,  3.8983930000e-05,
                1.6471520000e-04,  5.9859800000e-04,  1.9423900000e-03,
                5.6872370000e-03,  1.5013290000e-02,  3.4524550000e-02,
                6.4958200000e-02,  7.7161940000e-02,  -1.8734110000e-02,
                -3.0091850000e-01, -4.5546610000e-01, 1.2864630000e-01,
                7.1833160000e-01,  4.0517430000e-01,  2.1682270000e-02,
                -8.3435660000e-03, 3.6589790000e-03},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.0630080000e-07, -1.6041690000e-06, -8.4384370000e-06,
                -3.5631510000e-05, -1.2959980000e-04, -4.2015340000e-04,
                -1.2319540000e-03, -3.2489220000e-03, -7.4937170000e-03,
                -1.4101020000e-02, -1.6916000000e-02, 4.2189960000e-03,
                6.8338100000e-02,  1.0982010000e-01,  -4.0090050000e-02,
                -2.1747390000e-01, -2.4651350000e-01, 2.7314350000e-01,
                5.7483210000e-01,  3.0127130000e-01},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-4.0093670000e-07, -3.1892550000e-06, -1.6230790000e-05,
                -7.1579200000e-05, -2.4639580000e-04, -8.5449070000e-04,
                -2.3075930000e-03, -6.7282920000e-03, -1.3661650000e-02,
                -3.0622400000e-02, -2.6311370000e-02, -9.7601830000e-03,
                1.8019060000e-01,  1.5296340000e-01,  5.5054130000e-02,
                -9.5513640000e-01, 2.5868130000e-01,  1.8340490000e+00,
                -9.3332400000e-01, -6.9816050000e-01},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.9660420000e-07, -5.5680360000e-06, -2.8136840000e-05,
                -1.2524180000e-04, -4.2607870000e-04, -1.4990600000e-03,
                -3.9791030000e-03, -1.1856860000e-02, -2.3467340000e-02,
                -5.4677360000e-02, -4.3938200000e-02, -2.3761030000e-02,
                3.4359280000e-01,  3.1929600000e-01,  -1.3432070000e-01,
                -2.2210200000e+00, 2.5711420000e+00,  -2.2924040000e-01,
                -1.8324520000e+00, 1.5913330000e+00},
      doubles_t{4.3162650000e+06, 6.4634240000e+05, 1.4708970000e+05,
                4.1661520000e+04, 1.3590770000e+04, 4.9057500000e+03,
                1.9127460000e+03, 7.9260430000e+02, 3.4480650000e+02,
                1.5589990000e+02, 7.2230910000e+01, 3.2725060000e+01,
                1.5667620000e+01, 7.5034830000e+00, 3.3122230000e+00,
                1.5584710000e+00, 6.8391400000e-01, 1.4675700000e-01,
                7.0583000000e-02, 3.1449000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1449000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.6116000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2169000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.1000000000e-05, 3.6900000000e-04, 2.1290000000e-03,
                9.3690000000e-03, 3.3097000000e-02, 9.4431000000e-02,
                2.0807700000e-01, 3.3233300000e-01, 3.3298700000e-01,
                1.5684300000e-01, 2.1549000000e-02, -2.0950000000e-03,
                -1.7390000000e-03, -3.0000000000e-04, 2.9000000000e-05,
                -1.1000000000e-05},
      doubles_t{1.7745690000e+04, 4.2007210000e+03, 1.3644290000e+03,
                5.2208060000e+02, 2.2145950000e+02, 1.0090960000e+02,
                4.8401150000e+01, 2.3985360000e+01, 1.2182500000e+01,
                6.2422980000e+00, 3.1109440000e+00, 1.5099580000e+00,
                7.1084500000e-01, 2.7319000000e-01, 1.0423300000e-01,
                3.8291000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5000000000e-05, -1.3000000000e-04, -7.5100000000e-04,
                -3.3290000000e-03, -1.1912000000e-02, -3.4933000000e-02,
                -7.9989000000e-02, -1.3463600000e-01, -1.3859800000e-01,
                3.0278000000e-02, 3.3321600000e-01, 4.5615300000e-01,
                2.8505100000e-01, 4.6144000000e-02, -3.2490000000e-03,
                1.3570000000e-03},
      doubles_t{1.7745690000e+04, 4.2007210000e+03, 1.3644290000e+03,
                5.2208060000e+02, 2.2145950000e+02, 1.0090960000e+02,
                4.8401150000e+01, 2.3985360000e+01, 1.2182500000e+01,
                6.2422980000e+00, 3.1109440000e+00, 1.5099580000e+00,
                7.1084500000e-01, 2.7319000000e-01, 1.0423300000e-01,
                3.8291000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.1000000000e-05, 8.7000000000e-05, 5.4100000000e-04,
                2.2260000000e-03, 8.5930000000e-03, 2.3339000000e-02,
                5.8844000000e-02, 8.8289000000e-02, 1.2319200000e-01,
                -6.3186000000e-02, -3.5490200000e-01, -6.1970800000e-01,
                8.1298600000e-01, 8.1911800000e-01, -9.0170500000e-01,
                -1.3591300000e-01},
      doubles_t{1.7745690000e+04, 4.2007210000e+03, 1.3644290000e+03,
                5.2208060000e+02, 2.2145950000e+02, 1.0090960000e+02,
                4.8401150000e+01, 2.3985360000e+01, 1.2182500000e+01,
                6.2422980000e+00, 3.1109440000e+00, 1.5099580000e+00,
                7.1084500000e-01, 2.7319000000e-01, 1.0423300000e-01,
                3.8291000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{5.0000000000e-06, 4.2000000000e-05, 2.4100000000e-04,
                1.0850000000e-03, 3.8310000000e-03, 1.1423000000e-02,
                2.5792000000e-02, 4.4818000000e-02, 4.4598000000e-02,
                -1.1177000000e-02, -1.3813400000e-01, -1.8828500000e-01,
                -1.0739900000e-01, 4.4486300000e-01, 6.4023900000e-01,
                6.4457000000e-02},
      doubles_t{1.7745690000e+04, 4.2007210000e+03, 1.3644290000e+03,
                5.2208060000e+02, 2.2145950000e+02, 1.0090960000e+02,
                4.8401150000e+01, 2.3985360000e+01, 1.2182500000e+01,
                6.2422980000e+00, 3.1109440000e+00, 1.5099580000e+00,
                7.1084500000e-01, 2.7319000000e-01, 1.0423300000e-01,
                3.8291000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06, 2.9000000000e-05, 1.6500000000e-04,
                7.3400000000e-04, 2.6260000000e-03, 7.7250000000e-03,
                1.7733000000e-02, 3.0055000000e-02, 3.1094000000e-02,
                -1.0048000000e-02, -8.8306000000e-02, -1.2982400000e-01,
                -7.6937000000e-02, 2.1266100000e-01, 5.7306100000e-01,
                3.6965100000e-01},
      doubles_t{1.7745690000e+04, 4.2007210000e+03, 1.3644290000e+03,
                5.2208060000e+02, 2.2145950000e+02, 1.0090960000e+02,
                4.8401150000e+01, 2.3985360000e+01, 1.2182500000e+01,
                6.2422980000e+00, 3.1109440000e+00, 1.5099580000e+00,
                7.1084500000e-01, 2.7319000000e-01, 1.0423300000e-01,
                3.8291000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.8291000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{6.7110000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.6079000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.5326000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.1303000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.5300000000e-03, 2.5784000000e-02, 9.9119000000e-02,
                           2.3907300000e-01, 3.5719900000e-01, 3.6218800000e-01,
                           2.3646100000e-01, 6.0118000000e-02},
                 doubles_t{1.1334400000e+02, 3.3641400000e+01, 1.2331000000e+01,
                           4.9947800000e+00, 2.0728000000e+00, 8.3075300000e-01,
                           3.0917800000e-01, 1.0013000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.8900000000e-03, -2.8442000000e-02, -1.1242900000e-01,
                -2.7425700000e-01, -3.1554600000e-01, 5.7109000000e-02,
                5.6360400000e-01, 3.8463700000e-01},
      doubles_t{1.1334400000e+02, 3.3641400000e+01, 1.2331000000e+01,
                4.9947800000e+00, 2.0728000000e+00, 8.3075300000e-01,
                3.0917800000e-01, 1.0013000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{5.6950000000e-03, 4.2001000000e-02, 1.7354000000e-01,
                4.1015700000e-01, 1.1325200000e-01, -7.6968000000e-01,
                -3.1643000000e-02, 7.1379700000e-01},
      doubles_t{1.1334400000e+02, 3.3641400000e+01, 1.2331000000e+01,
                4.9947800000e+00, 2.0728000000e+00, 8.3075300000e-01,
                3.0917800000e-01, 1.0013000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0013000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{6.0522000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{2.2818000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{7.9200000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.5476000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0897000000e+00}));
    return abs_t(name, 26, r0, shells.begin(), shells.end());
} // cc_dash_pwcvtz_26

} // namespace chemcache
