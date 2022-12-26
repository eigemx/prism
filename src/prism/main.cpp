#include <fmt/color.h>
#include <fmt/core.h>
#include <unvpp/unvpp.h>


#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "mesh/from_unv.h"

auto main(int argc, char* argv[]) -> int {
    auto stream = std::ifstream("./external/unvpp/meshes/big_cube.unv");
    auto mesh = unv::read(stream);

    fmt::print("Count of vertices = {}\n", mesh.vertices.size());
    fmt::print("Count of elements = {}\n", mesh.elements.value().size());

    for (auto e : mesh.elements.value()) {
        if (e.vertices_ids.size() == 6) {
            auto faces = prism::mesh::hex_cell_faces(e.vertices_ids);

            for (auto face : faces) {
                for (auto i : face) {
                    fmt::print("{} ", i);
                }
                std::cout << std::endl;
            }
            break;
        }
    }
}