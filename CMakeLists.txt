cmake_minimum_required(VERSION 3.10)

project(
    prism
    LANGUAGES CXX
    VERSION 0.0.1
    DESCRIPTION "A Finite Volume Computational Fluid Dynamics Solver"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

add_executable(prism
    src/prism/main.cpp
    src/prism/types.h
    
    # mesh module
    src/prism/mesh/from_unv.cpp
    src/prism/mesh/face.cpp
    src/prism/mesh/cell.cpp
    
    src/prism/mesh/face.h
    src/prism/mesh/cell.h
    src/prism/mesh/from_unv.h

)

target_link_libraries(prism fmt unvpp)
