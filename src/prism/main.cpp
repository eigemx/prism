#include <fmt/color.h>
#include <fmt/core.h>
#include <toml++/toml.h>
#include <unvpp/unvpp.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string_view>
#include <tuple>
#include <utility>

#include "mesh/boundary.h"
#include "mesh/from_unv.h"

using namespace std::string_view_literals;

/*auto main(int argc, char* argv[]) -> int {
    std::string filename;

    if (argc >= 2) {
        filename = argv[1];
    } else {
        std::cout << "No input file found, using default test mesh\n";
        filename = "./external/unvpp/meshes/fine_cube.unv";
    }

    // test UNV library functionality
    auto pmesh = prism::mesh::UnvToPMesh(filename);
    pmesh.report_mesh_stats();

}*/

auto main(int argc, char* argv[]) -> int {
    std::string filename;

    if (argc >= 2) {
        filename = argv[1];
    } else {
        std::cout << "No input file found, using default test mesh\n";
        return -1;
    }
    auto bnames_to_look_for =
        std::vector<std::string_view> {"top", "bottom", "left", "right", "front", "back"};
    auto bcs = prism::mesh::read_boundary_conditions(filename, bnames_to_look_for);

    for (auto& bc : bcs) {
        std::cout << bc.name() << "\n";
    }
}
