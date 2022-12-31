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


auto main(int argc, char* argv[]) -> int {
    std::string filename;
    if (argc >= 2) {
        filename = argv[1];
    } else {
        std::cout << "No input file found, using default test mesh\n";
        filename = "./external/unvpp/meshes/fine_cube.unv";
    }
    // test UNV library functionality
    auto pmesh = prism::mesh::UnvToPMesh(filename);
    pmesh.quick_report();
}