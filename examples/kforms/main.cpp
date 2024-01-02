#include <fmt/core.h>
#include <prism/field.h>
#include <prism/mesh/unv.h>

#include <iterator>
#include <numeric>


auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    fmt::println("diffsolver - A steady state temperature diffusion solver");
    fmt::println("");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: diffsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "boundary.txt";
    fmt::print("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();
    fmt::println("Okay.");

    std::size_t bfaces_count {0};
    for (const auto& patch : mesh.boundary_patches()) {
        fmt::println("Patch `{}` has {} faces.", patch.name(), patch.faces_ids().size());
        bfaces_count += patch.faces_ids().size();
    }

    std::size_t correct_count =
        std::accumulate(std::begin(mesh.boundary_faces()),
                        std::end(mesh.boundary_faces()),
                        0,
                        [](int total, auto value) { return total + 1; } // NOLINT
        );

    fmt::println("Calculated count = {} | Correct count = {}", bfaces_count, correct_count);
}