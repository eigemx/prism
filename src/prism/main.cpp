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
    auto pmesh = prism::mesh::UnvToPMesh("./external/unvpp/meshes/fine_cube.unv");
}