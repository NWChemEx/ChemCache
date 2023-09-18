# Copyright 2022 NWChemEx-Project
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
set(VERSION 1.0.0) #TODO: Get from git
project(chemcache VERSION "${VERSION}" LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    nwx_cmake
    GIT_REPOSITORY https://github.com/NWChemEx-Project/NWXCMake
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(get_cmaize)
include(nwx_cxx_api_docs)

# Work out the project paths
set(project_inc_dir "${CMAKE_CURRENT_LIST_DIR}/include/${PROJECT_NAME}")
set(project_src_dir "${CMAKE_CURRENT_LIST_DIR}/src/${PROJECT_NAME}")

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)
option(BUILD_DOCS "Should we build the documentation?" OFF)
option(
    ONLY_BUILD_DOCS
    "If enabled and BUILD_DOCS is true no libraries will be built" OFF
)
option(BUILD_PYBIND11_PYBINDINGS "Use pybind11 to build Python3 bindings?" OFF)
option(BUILD_CPPYY_PYBINDINGS "Use Cppyy to build Python3 bindings?" OFF)

# Generate C++ API documentation, if BUILD_DOCS is set
nwx_cxx_api_docs("${project_inc_dir}")

cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx-Project/SimDE
    PRIVATE TRUE
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
)

cmaize_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${project_src_dir}"
    INCLUDE_DIRS "${project_inc_dir}"
    DEPENDS nwx::simde
)

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/python"
    DEPENDS "${PROJECT_NAME}"
)

include(nwx_python_mods)
cppyy_make_python_package(PACKAGE chemcache NAMESPACES chemcache DEPPACKAGES simde )

if("${BUILD_TESTING}")
    set(CXX_TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests/cxx")
    set(PYTHON_TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests/python")
    cmaize_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
	    VERSION v2.x
    )
    cmaize_add_tests(
        test_unit_${PROJECT_NAME}
        SOURCE_DIR "${CXX_TEST_DIR}/unit"
        INCLUDE_DIRS "${project_src_dir}"
        DEPENDS Catch2 ${PROJECT_NAME}
    )
    cmaize_add_tests(
        test_integration_${PROJECT_NAME}
        SOURCE_DIR "${CXX_TEST_DIR}/integration"
        INCLUDE_DIRS "${project_src_dir}"
        DEPENDS Catch2 ${PROJECT_NAME}
    )
    # Based on previous comment, not sure this test actually worked.
    #nwx_pybind11_tests(
    #scrape_bse ${PYTHON_TEST_DIR}/reference_data/test_scrape_bse.py
    #)
    nwx_pybind11_tests(
        py_chem_cache ${PYTHON_TEST_DIR}/unit_tests/test_chemcache.py
        SUBMODULES parallelzone pluginplay chemist simde
    )
endif()

cmaize_add_package(${PROJECT_NAME} NAMESPACE nwx::)
