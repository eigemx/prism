#include <prism/mesh/boundary.h>
#include <prism/mesh/unv.h>
#include <prism/print.h>
#include <toml++/toml.h>
#include <unvpp/unvpp.h>

#include <algorithm>
#include <cxxopts.hpp>
#include <filesystem>
#include <numeric>
#include <string>
#include <tabulate/tabulate.hpp>
#include <utility>
#include <vector>


void create_boundary_conditions_file(const std::filesystem::path& mesh_file) {
    // if mesh file extension is not .unv, warn user
    if (mesh_file.extension() != ".unv") {
        prism::warn(fmt::format("Input mesh file extension is not `.unv`."));
    }

    std::string bc_filename = "boundary.txt";

    // check if there are already a boundary conditions file in the same directory
    auto bc_filepath = mesh_file.parent_path() / bc_filename;

    if (std::filesystem::exists(bc_filepath)) {
        prism::error(
            fmt::format("Boundary conditions file `{}` already exists, please remove to proceed.",
                        bc_filename));
        return;
    }

    auto unv_mesh = unv::read(mesh_file);

    if (!unv_mesh.groups || unv_mesh.groups.value().empty()) {
        prism::error("No groups were found in mesh file");
        return;
    }

    toml::table table;

    for (const auto& group : unv_mesh.groups.value()) {
        table.insert(group.name(), toml::table {
                                       {"type", "wall"},
                                       {"velocity", toml::array {0.0, 0.0, 0.0}},
                                   });
    }

    std::ofstream file(bc_filepath);
    file << table;
    file.close();
}

void check_pmesh(const std::filesystem::path& mesh_file) {
    try {
        fmt::print("Loading mesh file: ");
        fmt::print(fg(fmt::color::dark_cyan), "{}\n\n", mesh_file.string());

        auto unv_mesh = prism::mesh::UnvToPMesh(mesh_file);

        // convert unv mesh to prism mesh
        auto prism_mesh = unv_mesh.to_pmesh();

        auto total_vol =
            std::accumulate(prism_mesh.cells().begin(), prism_mesh.cells().end(), 0.0,
                            [](double acc, const auto& cell) { return acc + cell.volume(); });

        auto surface_area = std::accumulate(prism_mesh.faces().begin(), prism_mesh.faces().end(),
                                            0.0, [](double acc, const auto& face) {
                                                if (face.has_neighbor()) {
                                                    return acc;
                                                }
                                                return acc + face.area();
                                            });

        fmt::print("Total volume  = {:.4f} L^3\n", total_vol);
        fmt::print("Surface area  = {:.4f} L^2\n", surface_area);


    } catch (const std::exception& e) {
        prism::error(e.what());
        return;
    }
}

auto main(int argc, char* argv[]) -> int {
    prism::print_header();

    // TODO: running with input file without argument does nothing, should print help

    try {
        cxxopts::Options options("pmesh", "A mesh utility for PMesh files");

        options.set_width(100).set_tab_expansion().add_options()
            // check input mesh and print quality stats
            ("c,check", "Check mesh file", cxxopts::value<std::string>())
            // create new boundary conditions file for the given input mesh, if not found
            ("n,new-boundary", "Create dummy boundary conditions file for the given input mesh",
             cxxopts::value<std::string>())
            // print help
            ("h,help", "Print help");

        auto result = options.parse(argc, argv);

        if (argc == 1) {
            fmt::print("Too few arguments, run with --help to see available options\n");
            return 0;
        }

        if (result.count("help") > 0) {
            fmt::print("{}\n", options.help());
            return 0;
        }

        if (result.count("check") > 0) {
            auto filename = result["check"].as<std::string>();
            check_pmesh(filename);
        }

        if (result.count("new-boundary") > 0) {
            auto filename = result["new-boundary"].as<std::string>();
            create_boundary_conditions_file(filename);
        }
    }

    catch (const cxxopts::exceptions::exception& e) {
        prism::error(fmt::format("Error parsing options: {}\n", e.what()));
        fmt::print("Run with --help to see available options\n");

        return 1;
    }
}
