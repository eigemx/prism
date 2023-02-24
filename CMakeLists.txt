cmake_minimum_required(VERSION 3.10)

project(
    prism
    LANGUAGES CXX
    VERSION 0.0.1
    DESCRIPTION "A Finite Volume Based Navier-Stokes Solver"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

# generate compile_commands.json to be used by clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

# setup tabulate
set(TABULATE_SOURCE_DIR "external/tabulate")
include_directories(${TABULATE_SOURCE_DIR}/single_include)

add_subdirectory(src/prism)
include_directories(src/)

# apps
add_subdirectory(src/apps/pmesh)