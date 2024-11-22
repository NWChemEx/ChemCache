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

abs_t cc_dash_pwcvtz_28() {
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
      doubles_t{8.2089960000e-06,  6.3828840000e-05,  3.3558000000e-04,
                1.4150750000e-03,  5.1244440000e-03,  1.6432560000e-02,
                4.6893980000e-02,  1.1635340000e-01,  2.3505110000e-01,
                3.3501840000e-01,  2.5347790000e-01,  7.3009010000e-02,
                6.1842440000e-02,  6.3029560000e-02,  1.0080630000e-02,
                -2.2445280000e-04, -5.9327670000e-05, -1.1585620000e-05,
                8.1151090000e-06,  -1.6816990000e-06},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.6578490000e-06, -2.8440940000e-05, -1.4959280000e-04,
                -6.3130090000e-04, -2.2930520000e-03, -7.4051230000e-03,
                -2.1520320000e-02, -5.5609740000e-02, -1.2301760000e-01,
                -2.1301040000e-01, -2.2658370000e-01, 3.5467960000e-02,
                5.1816970000e-01,  5.0256300000e-01,  8.9556740000e-02,
                -7.0313110000e-03, -4.3391670000e-04, -5.8317110000e-04,
                4.2287880000e-04,  -1.2667140000e-04},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{9.5941490000e-07,  7.4626140000e-06,  3.9238430000e-05,
                1.6578680000e-04,  6.0259050000e-04,  1.9556620000e-03,
                5.7303910000e-03,  1.5147560000e-02,  3.4934990000e-02,
                6.5980720000e-02,  7.8930830000e-02,  -1.9062490000e-02,
                -3.0959210000e-01, -4.5586100000e-01, 1.4829310000e-01,
                7.1340390000e-01,  3.9760630000e-01,  2.2955230000e-02,
                -9.1517580000e-03, 3.8754140000e-03},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.0137530000e-07, -1.5658320000e-06, -8.2371820000e-06,
                -3.4781050000e-05, -1.2652650000e-04, -4.1025890000e-04,
                -1.2038340000e-03, -3.1790620000e-03, -7.3538280000e-03,
                -1.3890220000e-02, -1.6778750000e-02, 4.1633780000e-03,
                6.8147030000e-02,  1.0610290000e-01,  -4.3399800000e-02,
                -2.0949500000e-01, -2.3102710000e-01, 2.5905320000e-01,
                5.6914260000e-01,  3.1581250000e-01},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-7.3036330000e-07, -5.8020130000e-06, -2.9585470000e-05,
                -1.3014990000e-04, -4.4953320000e-04, -1.5533140000e-03,
                -4.2198840000e-03, -1.2254190000e-02, -2.5197070000e-02,
                -5.6217170000e-02, -5.0222290000e-02, -1.6774120000e-02,
                3.3880210000e-01,  3.9849750000e-01,  -3.0320530000e-01,
                -2.0796190000e+00, 2.5005420000e+00,  -2.1690020000e-01,
                -1.7091780000e+00, 1.4691660000e+00},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.9242450000e-07, -3.1139090000e-06, -1.5904470000e-05,
                -6.9813940000e-05, -2.4178480000e-04, -8.3261950000e-04,
                -2.2702940000e-03, -6.5574270000e-03, -1.3542880000e-02,
                -2.9897680000e-02, -2.6931060000e-02, -7.8276930000e-03,
                1.7416670000e-01,  1.5954680000e-01,  1.9955500000e-02,
                -8.8970000000e-01, 2.4868920000e-01,  1.6130120000e+00,
                -5.9902770000e-01, -8.3690780000e-01},
      doubles_t{5.0459910000e+06, 7.5561420000e+05, 1.7195680000e+05,
                4.8704790000e+04, 1.5888410000e+04, 5.7351230000e+03,
                2.2361370000e+03, 9.2664680000e+02, 4.0317430000e+02,
                1.8234760000e+02, 8.4548850000e+01, 3.8396340000e+01,
                1.8458590000e+01, 8.8635480000e+00, 3.9162270000e+00,
                1.8388700000e+00, 8.0436200000e-01, 1.6979700000e-01,
                7.9306000000e-02, 3.4677000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{3.4677000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7408000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.5172000000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{4.1000000000e-05, 3.6300000000e-04, 2.0970000000e-03,
                9.2500000000e-03, 3.2796000000e-02, 9.4004000000e-02,
                2.0828000000e-01, 3.3365400000e-01, 3.3290400000e-01,
                1.5537200000e-01, 2.0859000000e-02, -2.4400000000e-03,
                -1.9980000000e-03, -3.3800000000e-04, 3.5000000000e-05,
                -1.2000000000e-05},
      doubles_t{2.1027920000e+04, 4.9775600000e+03, 1.6167400000e+03,
                6.1867180000e+02, 2.6251830000e+02, 1.1969070000e+02,
                5.7465850000e+01, 2.8528290000e+01, 1.4521480000e+01,
                7.4538500000e+00, 3.7235530000e+00, 1.8098130000e+00,
                8.5133600000e-01, 3.2481400000e-01, 1.1952200000e-01,
                4.2366000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.5000000000e-05, -1.2900000000e-04, -7.4900000000e-04,
                -3.3280000000e-03, -1.1947000000e-02, -3.5242000000e-02,
                -8.1204000000e-02, -1.3749300000e-01, -1.3922600000e-01,
                3.6016000000e-02, 3.3912800000e-01, 4.5047200000e-01,
                2.8178300000e-01, 4.7898000000e-02, -2.9870000000e-03,
                1.3090000000e-03},
      doubles_t{2.1027920000e+04, 4.9775600000e+03, 1.6167400000e+03,
                6.1867180000e+02, 2.6251830000e+02, 1.1969070000e+02,
                5.7465850000e+01, 2.8528290000e+01, 1.4521480000e+01,
                7.4538500000e+00, 3.7235530000e+00, 1.8098130000e+00,
                8.5133600000e-01, 3.2481400000e-01, 1.1952200000e-01,
                4.2366000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{6.0000000000e-06, 5.3000000000e-05, 3.0500000000e-04,
                1.3640000000e-03, 4.8760000000e-03, 1.4503000000e-02,
                3.3296000000e-02, 5.7482000000e-02, 5.8702000000e-02,
                -1.9904000000e-02, -1.9469500000e-01, -2.3961300000e-01,
                -2.2320000000e-03, 5.2143500000e-01, 5.4554000000e-01,
                4.3622000000e-02},
      doubles_t{2.1027920000e+04, 4.9775600000e+03, 1.6167400000e+03,
                6.1867180000e+02, 2.6251830000e+02, 1.1969070000e+02,
                5.7465850000e+01, 2.8528290000e+01, 1.4521480000e+01,
                7.4538500000e+00, 3.7235530000e+00, 1.8098130000e+00,
                8.5133600000e-01, 3.2481400000e-01, 1.1952200000e-01,
                4.2366000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.1000000000e-05, 9.5000000000e-05, 5.8000000000e-04,
                2.4510000000e-03, 9.2820000000e-03, 2.6009000000e-02,
                6.4096000000e-02, 1.0071000000e-01, 1.3253900000e-01,
                -6.5089000000e-02, -4.8975600000e-01, -4.9845500000e-01,
                9.6635700000e-01, 5.2837900000e-01, -8.6767600000e-01,
                -1.0445600000e-01},
      doubles_t{2.1027920000e+04, 4.9775600000e+03, 1.6167400000e+03,
                6.1867180000e+02, 2.6251830000e+02, 1.1969070000e+02,
                5.7465850000e+01, 2.8528290000e+01, 1.4521480000e+01,
                7.4538500000e+00, 3.7235530000e+00, 1.8098130000e+00,
                8.5133600000e-01, 3.2481400000e-01, 1.1952200000e-01,
                4.2366000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.0000000000e-06, 2.6000000000e-05, 1.5200000000e-04,
                6.7800000000e-04, 2.4270000000e-03, 7.2010000000e-03,
                1.6578000000e-02, 2.8392000000e-02, 2.8599000000e-02,
                -1.0132000000e-02, -8.2912000000e-02, -1.1599800000e-01,
                -7.2795000000e-02, 1.9564000000e-01, 5.6709900000e-01,
                3.9527000000e-01},
      doubles_t{2.1027920000e+04, 4.9775600000e+03, 1.6167400000e+03,
                6.1867180000e+02, 2.6251830000e+02, 1.1969070000e+02,
                5.7465850000e+01, 2.8528290000e+01, 1.4521480000e+01,
                7.4538500000e+00, 3.7235530000e+00, 1.8098130000e+00,
                8.5133600000e-01, 3.2481400000e-01, 1.1952200000e-01,
                4.2366000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.2366000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.6879000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{2.0080000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{5.7882000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{2.7127000000e+00}));
    shells.emplace_back(
      make_shell(pure_t::pure, 2,
                 doubles_t{3.3760000000e-03, 2.5141000000e-02, 9.7746000000e-02,
                           2.3470900000e-01, 3.4694500000e-01, 3.5106800000e-01,
                           2.5025500000e-01, 1.0008200000e-01},
                 doubles_t{1.4025270000e+02, 4.1726100000e+01, 1.5398100000e+01,
                           6.2771000000e+00, 2.6185000000e+00, 1.0526000000e+00,
                           3.9160000000e-01, 1.2620000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{-3.4950000000e-03, -2.6015000000e-02, -1.0387600000e-01,
                -2.5207000000e-01, -2.9458000000e-01, 1.1520000000e-03,
                4.3858900000e-01, 5.4362600000e-01},
      doubles_t{1.4025270000e+02, 4.1726100000e+01, 1.5398100000e+01,
                6.2771000000e+00, 2.6185000000e+00, 1.0526000000e+00,
                3.9160000000e-01, 1.2620000000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{5.0520000000e-03, 3.8087000000e-02, 1.5611300000e-01,
                3.8615600000e-01, 1.7560500000e-01, -6.2680700000e-01,
                -3.4277500000e-01, 7.9181900000e-01},
      doubles_t{1.4025270000e+02, 4.1726100000e+01, 1.5398100000e+01,
                6.2771000000e+00, 2.6185000000e+00, 1.0526000000e+00,
                3.9160000000e-01, 1.2620000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.2620000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{8.0748000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1196000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1286000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{7.3562000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.9771000000e+00}));
    return abs_t(name, 28, r0, shells.begin(), shells.end());
} // cc_dash_pwcvtz_28

} // namespace chemcache
