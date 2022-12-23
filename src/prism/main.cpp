#include <fmt/color.h>
#include <fmt/core.h>
#include <unvpp/unvpp.h>


#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>

auto main(int argc, char* argv[]) -> int {
    if (argc < 2) {
        fmt::print("Error: missing input mesh file\n");
        std::exit(-1);
    }

    auto stream = std::ifstream(argv[1]);
    auto mesh = unv::read(stream);

    fmt::print("Count of vertices = {}\n", mesh.vertices.size());
}