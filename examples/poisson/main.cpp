#include <fmt/core.h>
#include <prism/prism.h>

#include <string>
#include <vector>


auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    log::setLevel(log::Level::Debug);

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        log::error("Usage: poisson [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();

    auto P = field::Scalar("P", mesh, 0.0);


    return 0;
}
