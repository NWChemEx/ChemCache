# Copyright 2026 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)

# Downloads common CMake modules used throughout NWChemEx. Included before
# project() so nwx_set_version() is available to compute the version.
include(cmake/get_nwx_cmake.cmake)

# Version from scikit-build-core (wheels/sdists) or the latest git tag.
include(nwx_set_version)
nwx_set_version(chemcache_version "${CMAKE_CURRENT_LIST_DIR}")
project(chemcache VERSION "${chemcache_version}" LANGUAGES CXX)

include(disable_in_source_builds)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

option(BUILD_TESTING           "Whether to build the unit tests"    OFF)
option(BUILD_PYBIND11_BINDINGS "Build Python bindings via pybind11" ON)

# Documentation target (no-op unless BUILD_DOCS is on)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("cxx/include" "cxx/src")

### Dependencies ###
include(get_dependencies)
get_dependencies(simde)

### Library ###
include(nwx_library)
nwx_library(${PROJECT_NAME} "cxx/include" "cxx/src" simde)

### Python bindings (no-op unless BUILD_PYBIND11_BINDINGS is on) ###
include(nwx_python_module)
nwx_python_module(${PROJECT_NAME} "cxx/src")

### Tests ###
if(BUILD_TESTING)
    include(catch2_tests_from_dir)
    catch2_tests_from_dir(
        "test_unit_${PROJECT_NAME}"
        "tests/cxx/unit_tests"
        ${PROJECT_NAME}
        PRIVATE_INCLUDES "cxx/src"
    )

    include(nwx_python_test)
    nwx_python_test(
        py_unit_test_${PROJECT_NAME}
        "${CMAKE_CURRENT_LIST_DIR}/tests/python/unit_tests/test_chemcache.py"
    )
    nwx_python_test(
        py_utils_test_${PROJECT_NAME}
        "${CMAKE_CURRENT_LIST_DIR}/tests/python/utils_tests/test_utils.py"
    )
endif()
