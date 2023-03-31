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

#pragma once
#include <pluginplay/pluginplay.hpp>

namespace chemcache {

// Module declarations will go here
DECLARE_MODULE(elec_configs);

namespace elec_config_mods {

inline void set_defaults(pluginplay::ModuleManager& mm) {
    // Default submodules within this subcollection will be set here
}

inline void load_modules(pluginplay::ModuleManager& mm) {
    // Modules will be added to the ModuleManager here
    mm.add_module<elec_configs>("Elec Config from Z");

    set_defaults(mm);
}

} // namespace elec_config_mods

} // namespace chemcache
