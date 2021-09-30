cmake_minimum_required(VERSION 3.14)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION)
project(chemcache VERSION "${VERSION}" LANGUAGES CXX)

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(get_cpp)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}")

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)
option(BUILD_DOCS "Should we build the documentation?" OFF)
option(
    ONLY_BUILD_DOCS
    "If enabled and BUILD_DOCS is true no libraries will be built" OFF
)

option(BUILD_PYBINDINGS "Should we build Python3 bindings?" ON)

# Work out the project paths
set(project_inc_dir "${CMAKE_CURRENT_LIST_DIR}/include/${PROJECT_NAME}")
set(project_src_dir "${CMAKE_CURRENT_LIST_DIR}/src/${PROJECT_NAME}")

if("${BUILD_DOCS}")
    find_package(Doxygen REQUIRED)
    set(DOXYGEN_EXTRA_PACKAGES amsmath) # Extra LaTeX packages for documentation
    doxygen_add_docs(${PROJECT_NAME}_cxx_api "${project_inc_dir}")
    if("${ONLY_BUILD_DOCS}")
        return()
    endif()
endif()

cpp_find_or_build_dependency(
    libchemist
    URL github.com/NWChemEx-Project/libchemist
    PRIVATE TRUE
    BUILD_TARGET chemist
    FIND_TARGET nwx::chemist
    CMAKE_ARGS BUILD_TESTING=OFF
)

cpp_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${project_src_dir}"
    INCLUDE_DIR "${project_inc_dir}"
    DEPENDS nwx::chemist
)

if("${BUILD_TESTING}")
    set(project_test_dir "${CMAKE_CURRENT_LIST_DIR}/tests/${PROJECT_NAME}")
    cpp_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
	    VERSION v2.x
    )
    cpp_add_tests(
        test_${PROJECT_NAME}
        SOURCE_DIR "${project_test_dir}"
        INCLUDE_DIR "${project_src_dir}"
        DEPENDS Catch2::Catch2 ${PROJECT_NAME}
    )
    # Haven't gotten this to work yet.
    # add_test(NAME python_unittest 
    #          COMMAND "python -m unittest discover -s ${project_test_dir}")
endif()
