#include <fmt/core.h>
#include <prism/field.h>
#include <prism/mesh/unv.h>


auto main(int argc, char* argv[]) -> int {
    /*if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[0]); // NOLINT
        return -1;
    }

    auto mesh = prism::mesh::UnvToPMeshConverter(argv[1]); // NOLINT

    auto oneForm = prism::field::Vector()*/
    fmt::println("{} {}", sizeof(prism::Matrix3d), sizeof(double));
}