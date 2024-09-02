#TODO: review this file
cmake_minimum_required(VERSION 3.20)

project(
    prism
    LANGUAGES CXX
    VERSION 0.0.1
    DESCRIPTION "A C++ CFD library based on the finite volume method"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# check if prism is being the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(PRISM_MAIN_PROJECT ON)
else()
    set(PRISM_MAIN_PROJECT OFF)
endif()

find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)

# setup unvpp
set(UNVPP_SOURCE_DIR "external/unvpp")
add_subdirectory(${UNVPP_SOURCE_DIR})
include_directories(${UNVPP_SOURCE_DIR}/include)

# setup vtu11
set(VTU11_SOURCE_DIR "external/vtu11")
add_subdirectory(${VTU11_SOURCE_DIR})
include_directories(${VTU11_SOURCE_DIR}/)


# set build type to release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set PRISM_BUILD_EXAMPLES to ON by default
if(NOT DEFINED PRISM_BUILD_EXAMPLES)
    set(PRISM_BUILD_EXAMPLES ON)
    message(STATUS "PRISM_BUILD_EXAMPLES is set to ON by default, examples will be built.")
    message(STATUS "Set PRISM_BUILD_EXAMPLES to OFF to disable building examples.")
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
target_include_directories(prism PUBLIC src/)
target_link_libraries(prism 
                        PRIVATE Eigen3::Eigen fmt::fmt-header-only nlohmann_json::nlohmann_json unvpp
                        PUBLIC spdlog::spdlog)

if (PRISM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# TODO: tests should be used only when building as main project?
add_subdirectory(tests)
