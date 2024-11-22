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

#include "mini.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t mini_75() {
    // Basis Set name and origin point
    std::string name("mini");
    center_t r0(0.0, 0.0, 0.0);

    auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {
        return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
    };

    shells_t shells;
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{5.9972597752e-02, 3.6285678640e-01, 6.9492547395e-01},
      doubles_t{2.9263742000e+04, 4.4298015000e+03, 9.7511876000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-1.1237809258e-01, 7.3550235144e-01, 3.3133397812e-01},
      doubles_t{1.2932196000e+03, 1.2860003000e+02, 5.4807887000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-2.9625277406e-01, 9.0052002113e-01, 2.6977407637e-01},
      doubles_t{1.1845628000e+02, 2.3141641000e+01, 1.0379264000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.5506650695e-01, 8.4405001652e-01, 3.6499550714e-01},
      doubles_t{2.2971336000e+01, 5.3257089000e+00, 2.7498632000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{-3.1023300719e-01, 5.6258271303e-01, 6.0543011403e-01},
      doubles_t{4.7595836000e+00, 1.0178146000e+00, 5.0739080000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 0,
      doubles_t{1.0365820126e-01, 5.1251430625e-01, 4.8788210595e-01},
      doubles_t{9.1127990000e-01, 9.2000400000e-02, 3.4142700000e-02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{7.5559702640e-02, 4.1120171437e-01, 6.5129662276e-01},
      doubles_t{2.0244843000e+03, 4.7614657000e+02, 1.4333076000e+02}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{-3.7018801157e-02, 3.8857171214e-01, 6.7737642117e-01},
      doubles_t{2.7029118000e+02, 5.4486952000e+01, 2.2248123000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{3.6935547766e-01, 5.6372396590e-01, 1.3938349157e-01},
      doubles_t{1.0747778000e+01, 4.8975554000e+00, 2.2377669000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 1,
      doubles_t{1.9997509404e-01, 6.2258308143e-01, 2.6159259220e-01},
      doubles_t{1.8856196000e+00, 8.8924390000e-01, 3.4366700000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.0366219825e-01, 4.7100619207e-01, 5.9644248995e-01},
      doubles_t{3.1855059000e+02, 9.1648758000e+01, 3.0906581000e+01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{2.8105881645e-01, 5.5881543270e-01, 3.0997271814e-01},
      doubles_t{1.8298578000e+01, 7.1624476000e+00, 2.9489065000e+00}));
    shells.emplace_back(make_shell(
      pure_t::pure, 2,
      doubles_t{1.9619239896e-01, 5.4599669710e-01, 4.5014109761e-01},
      doubles_t{2.0604954000e+00, 7.4967700000e-01, 2.4248290000e-01}));
    shells.emplace_back(make_shell(
      pure_t::pure, 3,
      doubles_t{1.7838470901e-01, 5.3149592685e-01, 5.6523082855e-01},
      doubles_t{4.1648317000e+01, 1.2273287000e+01, 3.5704082000e+00}));
    return abs_t(name, 75, r0, shells.begin(), shells.end());
} // mini_75

} // namespace chemcache
