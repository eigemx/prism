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

#include "mesh/from_unv.h"

using namespace std::string_view_literals;

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
    pmesh.report_mesh_stats();

    static constexpr std::string_view some_toml = R"(
        [library]
        name = "toml++"
        authors = ["Mark Gillard <mark.gillard@outlook.com.au>"]
        cpp = 17
    )"sv;

    try {
        // parse directly from a string view:
        {
            toml::table tbl = toml::parse(some_toml);
            std::cout << tbl << "\n";
        }

    } catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }
}