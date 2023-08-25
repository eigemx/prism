cmake_minimum_required(VERSION 3.10)

project(
    prism
    LANGUAGES CXX
    VERSION 0.0.1
    DESCRIPTION "A Finite Volume Based CFD Framework"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# check if prism is being the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(PRISM_MAIN_PROJECT ON)
else()
    set(PRISM_MAIN_PROJECT OFF)
endif()

# setup eigen
set(EIGEN_SOURCE_DIR "external/eigen")
include_directories(${EIGEN_SOURCE_DIR})

# setup fmt
set(FMT_SOURCE_DIR "external/fmt")
add_subdirectory(${FMT_SOURCE_DIR})
include_directories(${FMT_SOURCE_DIR}/include)

# setup unvpp
set(UNVPP_SOURCE_DIR "external/unvpp")
add_subdirectory(${UNVPP_SOURCE_DIR})
include_directories(${UNVPP_SOURCE_DIR}/include)

# setup toml++
set(TOMLPP_SOURCE_DIR "external/tomlplusplus")
include_directories(${TOMLPP_SOURCE_DIR}/include)

# setup cxxopts
set(CXXOPTS_SOURCE_DIR "external/cxxopts")
include_directories(${CXXOPTS_SOURCE_DIR}/include)

# setup vtu11
set(VTU11_SOURCE_DIR "external/vtu11")
add_subdirectory(${VTU11_SOURCE_DIR})
include_directories(${VTU11_SOURCE_DIR}/)

# set build type to release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set SHOW_DEBUG_INFO to OFF by default
if(NOT DEFINED SHOW_DEBUG_INFO)
    set(SHOW_DEBUG_INFO OFF)
    message(STATUS "SHOW_DEBUG_INFO is not set, all debug info will be hidden.")
endif()

# if option SHOW_DEBUG_INFO is set, define PRISM_DEBUG
if(SHOW_DEBUG_INFO)
    message(STATUS "SHOW_DEBUG_INFO is set, all debug info will be shown.")
    add_definitions(-DPRISM_DEBUG)
endif()

# set optimization flags for release builds for gcc, clang and msvc
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    endif()
endif()

if (PRISM_MAIN_PROJECT)
    # generate compile_commands.json to be used by clang-tidy
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    # turn off asserts for release builds
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_definitions(-DNDEBUG)
    endif()

    # set output directories
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()


add_subdirectory(src/prism)
include_directories(src/)

# apps
add_subdirectory(src/apps/pmesh)

# solvers
add_subdirectory(src/apps/diffsolver)
add_subdirectory(src/apps/convsolver)