#include <fmt/color.h>
#include <fmt/core.h>
#include <unvpp/unvpp.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <utility>

#include "mesh/from_unv.h"


auto main() -> int {
    // test UNV library functionality
    auto mesh = unv::read("./external/unvpp/meshes/one_hex_cell.unv");

    fmt::print("Count of vertices = {}\n", mesh.vertices.size());
    fmt::print("Count of elements = {}\n", mesh.elements.value().size());
}