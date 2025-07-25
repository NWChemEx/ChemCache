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

#include "atoms/atoms.hpp"
#include "bases/bases.hpp"
#include "chemcache/chemcache_mm.hpp"
#include "density/density.hpp"
#include "hydrocarbons/hydrocarbons.hpp"
#include "molecules/molecules.hpp"

namespace chemcache {

inline void set_defaults(pluginplay::ModuleManager& mm) {
    // Default submodules between collections can be set here
    mm.change_submod("NWX Molecules", "Atoms", "Atom");
}

void load_modules(pluginplay::ModuleManager& mm) {
    // Add subcollection load calls here
    atom_mods::load_modules(mm);
    bases_mods::load_modules(mm);
    density_mods::load_modules(mm);
    hydrocarbons::load_modules(mm);
    molecule_mods::load_modules(mm);

    set_defaults(mm);
}

} // namespace chemcache
