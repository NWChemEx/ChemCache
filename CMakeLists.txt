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

#Downloads common CMake modules used throughout NWChemEx
include(cmake/get_nwx_cmake.cmake)

#Sets the version to whatever git thinks it is
include(get_version_from_git)
get_version_from_git(chemcache_version "${CMAKE_CURRENT_LIST_DIR}")
project(chemcache VERSION "${chemcache_version}" LANGUAGES CXX)

include(nwx_versions)
include(get_cmaize)
include(nwx_cxx_api_docs)

### Paths ###
set(CHEMCACHE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CHEMCACHE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(CHEMCACHE_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
set(CHEMCACHE_EXP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/experimental")

nwx_cxx_api_docs("${CHEMCACHE_SOURCE_DIR}" "${CHEMCACHE_INCLUDE_DIR}")

### Options ###
cmaize_option_list(
    BUILD_TESTING OFF "Should we build the tests?"
    BUILD_PYBIND11_PYBINDINGS ON "Build pybind11 python3 bindings?"
    ENABLE_EXPERIMENTAL_FEATURES OFF "Build features which are not 1.0-ed yet?"
)

cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx/SimDE
    PRIVATE TRUE
    VERSION ${NWX_SIMDE_VERSION}
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
               BUILD_PYBIND11_PYBINDINGS=${BUILD_PYBIND11_PYBINDINGS}
               ENABLE_EXPERIMENTAL_FEATURES=${ENABLE_EXPERIMENTAL_FEATURES}
)

cmaize_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${CHEMCACHE_SOURCE_DIR}/${PROJECT_NAME}"
    INCLUDE_DIRS "${CHEMCACHE_INCLUDE_DIR}/${PROJECT_NAME}"
    DEPENDS simde
)

if("${ENABLE_EXPERIMENTAL_FEATURES}")

    cmaize_add_library(
        experimental_${PROJECT_NAME}
        SOURCE_DIR "${CHEMCACHE_EXP_DIR}/src/${PROJECT_NAME}"
        INCLUDE_DIRS "${CHEMCACHE_EXP_DIR}/include/${PROJECT_NAME}"
        DEPENDS ${PROJECT_NAME}
    )

else()
    add_library(experimental_${PROJECT_NAME} INTERFACE)
endif()

# N.B. this is a no-op if BUILD_PYBIND11_PYBINDINGS is not turned on
include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${CHEMCACHE_SOURCE_DIR}/python"
    DEPENDS "${PROJECT_NAME}"
)

if("${BUILD_TESTING}")
    set(CXX_TEST_DIR ${CHEMCACHE_TESTS_DIR}/cxx)
    set(PYTHON_TEST_DIR ${CHEMCACHE_TESTS_DIR}/python)

    cmaize_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
	    VERSION ${NWX_CATCH2_VERSION}
    )

    cmaize_add_tests(
        test_unit_${PROJECT_NAME}
        SOURCE_DIR "${CXX_TEST_DIR}/unit_tests"
        INCLUDE_DIRS "${CHEMCACHE_SOURCE_DIR}/${PROJECT_NAME}"
        DEPENDS Catch2 ${PROJECT_NAME}
    )

    # N.B. these are no-ops if BUILD_PYBIND11_PYBINDINGS is not turned on
    nwx_pybind11_tests(
        py_chem_cache "${PYTHON_TEST_DIR}/unit_tests/test_chemcache.py"
        SUBMODULES parallelzone pluginplay chemist simde
    )

    if("${ENABLE_EXPERIMENTAL_FEATURES}")

        cmaize_add_tests(
            test_unit_experimental_${PROJECT_NAME}
            SOURCE_DIR "${CHEMCACHE_EXP_DIR}/tests/cxx/unit_tests"
            INCLUDE_DIRS "${CHEMCACHE_EXP_DIR}/src/${PROJECT_NAME}"
            DEPENDS Catch2::Catch2 experimental_${PROJECT_NAME}
        )

        nwx_pybind11_tests(
            py_experimental_${PROJECT_NAME}
            "${CHEMCACHE_EXP_DIR}/tests/python/unit_tests/test_python.py"
            SUBMODULES pluginplay chemist parallelzone
        )

    endif()

endif()

cmaize_add_package(
    ${PROJECT_NAME} experimental_${PROJECT_NAME} NAMESPACE nwx::
)
