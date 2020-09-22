cmake_minimum_required(VERSION 3.14)
project(ReferenceData VERSION 1.0.0)

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(get_cpp)

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)

cpp_find_or_build_dependency(
    chemist
    URL github.com/NWChemEx-Project/libchemist
    PRIVATE TRUE
    BUILD_TARGET chemist
    FIND_TARGET nwx::chemist
    CMAKE_ARGS BUILD_TESTING=OFF
)
cpp_add_library(
    reference_data
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/libchemist"
    INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/libchemist"
    DEPENDS nwx::chemist
)

if("${BUILD_TESTING}")
    cpp_add_tests(
        test_reference_data
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/tests"
        INCLUDE_DIR ""
        DEPENDS reference_data
    )
endif()


