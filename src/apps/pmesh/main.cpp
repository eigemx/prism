#include <prism/core.h>
#include <toml++/toml.h>
#include <unvpp/unvpp.h>

#include <algorithm>
#include <cxxopts.hpp>
#include <filesystem>
#include <numeric>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "export_vtu.h"

void create_boundary_conditions_file(const std::filesystem::path& mesh_file) {
    // if mesh file extension is not .unv, warn user
    if (mesh_file.extension() != ".unv") {
        prism::warn(prism::format("Input mesh file extension is not `.unv`."));
    }

    std::string bc_filename = "boundary.txt";

    // check if there are already a boundary conditions file in the same directory
    auto bc_filepath = mesh_file.parent_path() / bc_filename;

    if (std::filesystem::exists(bc_filepath)) {
        prism::error(prism::format(
            "Boundary conditions file `{}` already exists, please remove to proceed.",
            bc_filename));
        return;
    }

    auto unv_mesh = unvpp::read(mesh_file);

    if (!unv_mesh.groups || unv_mesh.groups.value().empty()) {
        prism::error("No groups were found in mesh file");
        return;
    }

    prism::print("Creating boundary conditions file...\n");

    toml::table table;

    for (const auto& group : unv_mesh.groups.value()) {
        table.insert(group.name(),
                     toml::table {
                         {"type", "wall"},
                         {"phi", 0.0},
                     });
    }

    std::ofstream file(bc_filepath);
    file << table;
    file.close();
}

auto min_max_face_area(const prism::mesh::PMesh& prism_mesh) {
    auto min_max = std::minmax_element(
        prism_mesh.faces().begin(),
        prism_mesh.faces().end(),
        [](const auto& face1, const auto& face2) { return face1.area() < face2.area(); });

    return std::make_pair(min_max.first->area(), min_max.second->area());
}

auto min_max_face_non_ortho(const prism::mesh::PMesh& prism_mesh) {
    // return min and max non orthogonality of faces for internal faces only
    auto min_non_ortho = std::numeric_limits<double>::max();
    auto max_non_ortho = std::numeric_limits<double>::min();

    for (const auto& face : prism_mesh.faces()) {
        if (!face.has_neighbor()) {
            continue;
        }

        auto non_ortho = prism_mesh.face_non_ortho(face);

        if (non_ortho < min_non_ortho) {
            min_non_ortho = non_ortho;
        }

        if (non_ortho > max_non_ortho) {
            max_non_ortho = non_ortho;
        }
    }

    return std::make_pair(min_non_ortho, max_non_ortho);
}


void check_pmesh(const std::filesystem::path& mesh_file) {
    try {
        prism::print("Loading mesh file: ");
        prism::print(fg(fmt::color::dark_cyan), "`{}`...\n\n", mesh_file.string());

        auto unv_mesh = prism::mesh::UnvToPMesh(mesh_file);

        // convert unv mesh to prism mesh
        auto prism_mesh = unv_mesh.to_pmesh();

        auto n_boundary_faces = std::count_if(
            prism_mesh.faces().begin(), prism_mesh.faces().end(), [](const auto& face) {
                return !face.has_neighbor();
            });

        prism::print("Mesh Elements:\n");
        prism::print(" - Vertices count = {} vertices\n", prism_mesh.vertices().size());
        prism::print(" - Faces count =  {} faces (including {} boundary faces)\n",
                     prism_mesh.faces().size(),
                     n_boundary_faces);
        prism::print(" - Cells count = {} cells\n\n", prism_mesh.cells().size());

        auto total_vol =
            std::accumulate(prism_mesh.cells().begin(),
                            prism_mesh.cells().end(),
                            0.0,
                            [](double acc, const auto& cell) { return acc + cell.volume(); });

        auto surface_area = std::accumulate(prism_mesh.faces().begin(),
                                            prism_mesh.faces().end(),
                                            0.0,
                                            [](double acc, const auto& face) {
                                                if (face.has_neighbor()) {
                                                    return acc;
                                                }
                                                return acc + face.area();
                                            });

        prism::print("Total volume  = {:.4f} L^3\n", total_vol);
        prism::print("Surface area  = {:.4f} L^2\n\n", surface_area);

        auto [min_face_area, max_face_area] = min_max_face_area(prism_mesh);
        auto [min_face_non_ortho, max_face_non_ortho] = min_max_face_non_ortho(prism_mesh);

        prism::print("Max face area = {:.4f} L^2 and min face area = {:.4f} L^2\n",
                     max_face_area,
                     min_face_area);

        prism::print(
            "Max face non orthogonality = {:.4f} and min face non orthogonality = {:.4f}\n",
            max_face_non_ortho,
            min_face_non_ortho);

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
            ("c,check", "Check mesh file", cxxopts::value<std::string>())(
                "e,export-vtu", "Export pmesh to vtu file", cxxopts::value<std::string>())
            // create new boundary conditions file for the given input mesh, if not found
            ("n,new-boundary",
             "Create dummy boundary conditions file for the given input mesh",
             cxxopts::value<std::string>())
            // print help
            ("h,help", "Print help");


        auto result = options.parse(argc, argv);

        if (argc == 1) {
            prism::print("Too few arguments, run with --help to see available options\n");
            return 0;
        }

        if (result.count("help") > 0) {
            prism::print("{}\n", options.help());
            return 0;
        }

        if (result.count("check") > 0) {
            auto filename = result["check"].as<std::string>();
            check_pmesh(filename);
        }

        if (result.count("new-boundary") > 0) {
            auto filename = result["new-boundary"].as<std::string>();
            create_boundary_conditions_file(filename);
            prism::print("File `boundary.txt` was created successfully for mesh file: {}\n",
                         filename);
        }

        if (result.count("export-vtu") > 0) {
            // TODO: clean this mess
            auto filename = result["export-vtu"].as<std::string>();

            // remove .unv extension if exists and add .vtu extension
            filename = std::regex_replace(filename, std::regex("\\.unv$"), "");
            filename += ".vtu";

            prism::print("Loading mesh file: ");
            prism::print(
                fg(fmt::color::dark_cyan), "`{}`...\n", result["export-vtu"].as<std::string>());

            auto unv_mesh = prism::mesh::UnvToPMesh(result["export-vtu"].as<std::string>());

            // convert unv mesh to prism mesh
            auto prism_mesh = unv_mesh.to_pmesh();

            prism::print("Exporting mesh to vtu file...\n");

            export_to_vtu(prism_mesh, filename);
            prism::print("File `{}` was created successfully\n", filename);
        }
    }

    catch (const cxxopts::exceptions::exception& e) {
        prism::error(prism::format("Error parsing options: {}\n", e.what()));
        prism::print("Run with --help to see available options\n");

        return 1;
    }
}
