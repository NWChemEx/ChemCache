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

PLUGIN_CTOR(ChemCache) {
    set_name("chemcache");

    // The defaults refer to modules from more than one subcollection, so all
    // modules are added before any defaults are set
    atom_mods::load_modules(modules());
    bases_mods::load_modules(modules());
    density_mods::load_modules(modules());
    hydrocarbons::load_modules(modules());
    molecule_mods::load_modules(modules());

    atom_mods::set_defaults(modules());
    bases_mods::set_defaults(modules());
    density_mods::set_defaults(modules());
    hydrocarbons::set_defaults(modules());
    molecule_mods::set_defaults(modules());
    change_submod("NWX Molecules", "Atoms", "Atom");
}

} // namespace chemcache
