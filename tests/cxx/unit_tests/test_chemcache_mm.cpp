/*
 * Copyright 2026 NWChemEx-Project
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

#include "catch.hpp"
#include "chemcache/chemcache.hpp"

TEST_CASE("ChemCache plugin") {
    chemcache::ChemCache plugin;
    const auto& modules = plugin.get_modules();

    SECTION("Identity") {
        REQUIRE(plugin.name() == "chemcache");
        // Every default comes from the plugin itself
        REQUIRE(plugin.preferred_providers().empty());
    }

    SECTION("Modules from each subcollection") {
        REQUIRE(modules.count("Atom"));
        REQUIRE(modules.count("sto-3g"));
        REQUIRE(modules.count("sto-3g SAD density"));
        REQUIRE(modules.count("Straight Chain Alkane"));
        REQUIRE(modules.count("NWX Molecules"));
    }

    SECTION("Defaults set within and across subcollections") {
        REQUIRE(modules.has_submod("sto-3g", "Atomic Basis"));
        REQUIRE(modules.has_submod("sto-3g SAD density", "Atomic Basis"));
        REQUIRE(modules.has_submod("sto-3g SAD density", "Atomic Density"));
        REQUIRE(modules.has_submod("NWX Molecules", "Atoms"));
    }

    SECTION("Can be loaded into a ModuleManager") {
        pluginplay::ModuleManager mm;
        mm.load_plugin<chemcache::ChemCache>();
        REQUIRE(mm.has_plugin("chemcache"));
        REQUIRE(mm.size() == modules.size());
        REQUIRE(mm.at("NWX Molecules").submods().at("Atoms").value() ==
                mm.at("Atom"));
    }
}
