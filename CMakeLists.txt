cmake_minimum_required(VERSION 3.20)

project(
    prism
    LANGUAGES CXX
    VERSION 0.0.1
    DESCRIPTION "A C++ CFD library based on the finite volume method"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(PRISM_MAIN_PROJECT ON)
else()
    set(PRISM_MAIN_PROJECT OFF)
endif()

# Dependencies
find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)
add_subdirectory(external/unvpp)
add_subdirectory(external/vtu11)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "Setting build type to 'Release' as none was specified.")
endif()

# set optimization flags for release builds for gcc, clang and msvc
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    endif()
endif()

# Main project specific setup
if (PRISM_MAIN_PROJECT)
    # Generate compile_commands.json for tooling (e.g. clangd)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    # Set output directories for executables and libraries
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

    # Option to build examples, ON by default
    option(PRISM_BUILD_EXAMPLES "Build the examples" ON)
    if (PRISM_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()

    # Add tests only when building as the main project
    add_subdirectory(tests)
endif()

# Main library target
add_subdirectory(src/prism)

target_include_directories(prism
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

target_link_libraries(prism
    PUBLIC
        Eigen3::Eigen
        fmt::fmt-header-only
        nlohmann_json::nlohmann_json
        unvpp
        vtu11
        spdlog::spdlog
)
