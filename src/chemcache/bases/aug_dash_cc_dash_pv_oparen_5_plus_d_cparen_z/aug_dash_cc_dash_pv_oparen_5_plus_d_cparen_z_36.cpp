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

#include "aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z_36() {
    // Basis Set name and origin point
    std::string name("aug-cc-pv(5+d)z");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.8000000000e-07,  1.4200000000e-06,  7.4600000000e-06,
                3.1590000000e-05,  1.1575000000e-04,  3.8150000000e-04,
                1.1590000000e-03,  3.2934000000e-03,  8.8161000000e-03,
                2.2218000000e-02,  5.2088100000e-02,  1.1063560000e-01,
                2.0253260000e-01,  2.9263500000e-01,  2.8512240000e-01,
                1.4550640000e-01,  2.3993900000e-02,  -9.4900000000e-05,
                5.5780000000e-04,  -2.8700000000e-04, 9.6600000000e-05,
                -7.8400000000e-05, 3.8770000000e-05,  -1.8220000000e-05,
                7.6300000000e-06,  -1.5900000000e-06},
      doubles_t{
        1.8282210000e+08, 2.7356160000e+07, 6.2211700000e+06, 1.7602780000e+06,
        5.7319380000e+05, 2.0625850000e+05, 8.0026670000e+04, 3.2939080000e+04,
        1.4222630000e+04, 6.3930710000e+03, 2.9764540000e+03, 1.4305250000e+03,
        7.0792620000e+02, 3.5984850000e+02, 1.8714970000e+02, 9.8634520000e+01,
        5.0547870000e+01, 2.7167000000e+01, 1.4615100000e+01, 7.6513520000e+00,
        3.9972630000e+00, 2.0858530000e+00, 1.0147970000e+00, 5.1978800000e-01,
        2.4510300000e-01, 1.1189600000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-5.7000000000e-07, -4.5000000000e-07, -2.3500000000e-06,
                -9.9400000000e-06, -3.6400000000e-05, -1.2010000000e-04,
                -3.6540000000e-04, -1.0407000000e-03, -2.8038000000e-03,
                -7.1509000000e-03, -1.7220400000e-02, -3.8480000000e-02,
                -7.7862800000e-02, -1.3474230000e-01, -1.7761480000e-01,
                -1.0684130000e-01, 1.8961320000e-01,  5.0918710000e-01,
                3.9398590000e-01,  9.1903200000e-02,  3.9195000000e-03,
                1.7496000000e-03,  -4.2580000000e-04, 1.5290000000e-04,
                -9.3810000000e-05, 1.3260000000e-05},
      doubles_t{
        1.8282210000e+08, 2.7356160000e+07, 6.2211700000e+06, 1.7602780000e+06,
        5.7319380000e+05, 2.0625850000e+05, 8.0026670000e+04, 3.2939080000e+04,
        1.4222630000e+04, 6.3930710000e+03, 2.9764540000e+03, 1.4305250000e+03,
        7.0792620000e+02, 3.5984850000e+02, 1.8714970000e+02, 9.8634520000e+01,
        5.0547870000e+01, 2.7167000000e+01, 1.4615100000e+01, 7.6513520000e+00,
        3.9972630000e+00, 2.0858530000e+00, 1.0147970000e+00, 5.1978800000e-01,
        2.4510300000e-01, 1.1189600000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{2.0000000000e-08,  1.8000000000e-07,  9.4000000000e-07,
                3.9900000000e-06,  1.4620000000e-05,  4.8190000000e-05,
                1.4670000000e-04,  4.1770000000e-04,  1.1267000000e-03,
                2.8769000000e-03,  6.9549000000e-03,  1.5636500000e-02,
                3.2080100000e-02,  5.6868900000e-02,  7.8484500000e-02,
                5.0339800000e-02,  -1.0427420000e-01, -3.7437610000e-01,
                -4.0111310000e-01, 9.6838800000e-02,  6.4287760000e-01,
                4.8606000000e-01,  7.3345600000e-02,  -1.0077000000e-03,
                1.6087000000e-03,  -4.3930000000e-05},
      doubles_t{
        1.8282210000e+08, 2.7356160000e+07, 6.2211700000e+06, 1.7602780000e+06,
        5.7319380000e+05, 2.0625850000e+05, 8.0026670000e+04, 3.2939080000e+04,
        1.4222630000e+04, 6.3930710000e+03, 2.9764540000e+03, 1.4305250000e+03,
        7.0792620000e+02, 3.5984850000e+02, 1.8714970000e+02, 9.8634520000e+01,
        5.0547870000e+01, 2.7167000000e+01, 1.4615100000e+01, 7.6513520000e+00,
        3.9972630000e+00, 2.0858530000e+00, 1.0147970000e+00, 5.1978800000e-01,
        2.4510300000e-01, 1.1189600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0147970000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{5.1978800000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-6.0000000000e-09, -5.7000000000e-08, -3.0000000000e-07,
                -1.2700000000e-06, -4.6580000000e-06, -1.5347000000e-05,
                -4.6736000000e-05, -1.3302200000e-04, -3.5906200000e-04,
                -9.1650000000e-04, -2.2184000000e-03, -4.9883000000e-03,
                -1.0266100000e-02, -1.8244500000e-02, -2.5411000000e-02,
                -1.6393100000e-02, 3.4697700000e-02,  1.3212830000e-01,
                1.4709250000e-01,  -4.1821600000e-02, -3.1952400000e-01,
                -4.3632860000e-01, 4.2911000000e-03,  5.6455580000e-01,
                5.5028450000e-01,  1.2977100000e-01},
      doubles_t{
        1.8282210000e+08, 2.7356160000e+07, 6.2211700000e+06, 1.7602780000e+06,
        5.7319380000e+05, 2.0625850000e+05, 8.0026670000e+04, 3.2939080000e+04,
        1.4222630000e+04, 6.3930710000e+03, 2.9764540000e+03, 1.4305250000e+03,
        7.0792620000e+02, 3.5984850000e+02, 1.8714970000e+02, 9.8634520000e+01,
        5.0547870000e+01, 2.7167000000e+01, 1.4615100000e+01, 7.6513520000e+00,
        3.9972630000e+00, 2.0858530000e+00, 1.0147970000e+00, 5.1978800000e-01,
        2.4510300000e-01, 1.1189600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{2.4510300000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1189600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 0, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4277000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{2.9700000000e-05, 2.6510000000e-04, 1.5416000000e-03,
                6.9065000000e-03, 2.5139700000e-02, 7.5012400000e-02,
                1.7674330000e-01, 3.0751350000e-01, 3.4706440000e-01,
                2.0028020000e-01, 4.3050800000e-02, 2.4772000000e-03,
                6.7890000000e-04, -8.2100000000e-05, 4.6900000000e-05,
                -2.2100000000e-05, 5.8000000000e-06},
      doubles_t{4.2993060000e+04, 1.0173720000e+04, 3.3031060000e+03,
                1.2635400000e+03, 5.3636550000e+02, 2.4487620000e+02,
                1.1799120000e+02, 5.9021250000e+01, 3.0356070000e+01,
                1.5819980000e+01, 8.1045800000e+00, 4.0979640000e+00,
                2.0560610000e+00, 9.5214500000e-01, 4.4477400000e-01,
                1.9749600000e-01, 8.3823000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-1.2100000000e-05, -1.0780000000e-04, -6.2900000000e-04,
                -2.8323000000e-03, -1.0446200000e-02, -3.1940000000e-02,
                -7.8459900000e-02, -1.4397190000e-01, -1.6917030000e-01,
                -1.7596600000e-02, 3.0026490000e-01, 4.8476610000e-01,
                2.9672480000e-01, 4.7453300000e-02, -4.7280000000e-04,
                9.2530000000e-04, -1.8470000000e-04},
      doubles_t{4.2993060000e+04, 1.0173720000e+04, 3.3031060000e+03,
                1.2635400000e+03, 5.3636550000e+02, 2.4487620000e+02,
                1.1799120000e+02, 5.9021250000e+01, 3.0356070000e+01,
                1.5819980000e+01, 8.1045800000e+00, 4.0979640000e+00,
                2.0560610000e+00, 9.5214500000e-01, 4.4477400000e-01,
                1.9749600000e-01, 8.3823000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5214500000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{4.4477400000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.3000000000e-06, 2.9300000000e-05, 1.7130000000e-04,
                7.6950000000e-04, 2.8514000000e-03, 8.7204000000e-03,
                2.1618100000e-02, 3.9802400000e-02, 4.7477500000e-02,
                -4.7730000000e-04, -1.0218910000e-01, -1.8236110000e-01,
                -1.1733630000e-01, 1.9542440000e-01, 4.5901760000e-01,
                4.0961740000e-01, 1.2597150000e-01},
      doubles_t{4.2993060000e+04, 1.0173720000e+04, 3.3031060000e+03,
                1.2635400000e+03, 5.3636550000e+02, 2.4487620000e+02,
                1.1799120000e+02, 5.9021250000e+01, 3.0356070000e+01,
                1.5819980000e+01, 8.1045800000e+00, 4.0979640000e+00,
                2.0560610000e+00, 9.5214500000e-01, 4.4477400000e-01,
                1.9749600000e-01, 8.3823000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9749600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{8.3823000000e-02}));
    shells.emplace_back(make_shell(pure_t::pure, 1, doubles_t{1.0000000000e+00},
                                   doubles_t{3.3129000000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{4.9600000000e-05, 4.9440000000e-04, 3.0265000000e-03,
                1.3346100000e-02, 4.3786900000e-02, 1.1143880000e-01,
                2.1303410000e-01, 2.9792410000e-01, 3.0796600000e-01,
                2.1107750000e-01, 7.6326100000e-02, 9.1733000000e-03},
      doubles_t{2.0674360000e+03, 6.2569370000e+02, 2.4394680000e+02,
                1.0842370000e+02, 5.2005220000e+01, 2.6115400000e+01,
                1.3546750000e+01, 7.1058100000e+00, 3.7215540000e+00,
                1.9291200000e+00, 9.5582600000e-01, 4.0519700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.9291200000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{9.5582600000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{4.0519700000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7410000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 2, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0140000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.0940000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8700000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{3.1500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 3, doubles_t{1.0000000000e+00},
                                   doubles_t{1.7840000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{1.1040000000e+00}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{5.0100000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 4, doubles_t{1.0000000000e+00},
                                   doubles_t{2.5500000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{9.3030000000e-01}));
    shells.emplace_back(make_shell(pure_t::pure, 5, doubles_t{1.0000000000e+00},
                                   doubles_t{5.8000000000e-01}));
    return abs_t(name, 36, r0, shells.begin(), shells.end());
} // aug_dash_cc_dash_pv_oparen_5_plus_d_cparen_z_36

} // namespace chemcache
