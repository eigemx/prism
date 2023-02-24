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
#include <utility>
#include <vector>


auto boundary_type_to_string(prism::mesh::BoundaryConditionType type) -> std::string {
    switch (type) {
        case prism::mesh::BoundaryConditionType::Wall:
            return "wall";
        case prism::mesh::BoundaryConditionType::Inlet:
            return "inlet";
        case prism::mesh::BoundaryConditionType::Outlet:
            return "outlet";
        case prism::mesh::BoundaryConditionType::Gradient:
            return "gradient";
        case prism::mesh::BoundaryConditionType::Symmetry:
            return "symmetry";
        case prism::mesh::BoundaryConditionType::Empty:
            return "empty";
        case prism::mesh::BoundaryConditionType::Unknown:
            return "unknown";
    }
}

template <typename T>
auto min_max(const std::vector<T>& vec) -> std::pair<T, T> {
    auto min = std::numeric_limits<T>::infinity();
    auto max = -std::numeric_limits<T>::infinity();

    for (const auto& val : vec) {
        if (val < min) {
            min = val;
        }

        if (val > max) {
            max = val;
        }
    }

    return {min, max};
}

void create_boundary_conditions_file(const std::filesystem::path& mesh_file) {
    std::string filename = "boundary.txt";

    // check if there are already a boundary conditions file in the same directory
    auto boundary_file = mesh_file.parent_path() / filename;

    if (std::filesystem::exists(boundary_file)) {
        prism::error(fmt::format(
            "Boundary conditions file `{}` already exists, please remove to proceed.", filename));
        return;
    }

    auto unv_mesh = unv::read(mesh_file);

    if (!unv_mesh.groups || unv_mesh.groups.value().empty()) {
        prism::error("No groups found in mesh file");
        return;
    }

    toml::table table;

    for (const auto& group : unv_mesh.groups.value()) {
        table.insert(group.name(), toml::table {
                                       {"type", "wall"},
                                       {"velocity", toml::array {0.0, 0.0, 0.0}},
                                   });
    }

    std::ofstream file(boundary_file);
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

        // print number of cells and faces
        fmt::print("Elements stats:\n");
        fmt::print("---------------\n");
        fmt::print("Number of cells = {} cell\n", prism_mesh.cells().size());
        fmt::print("Number of faces = {} face\n\n", prism_mesh.faces().size());

        // calculate total volume of mesh
        auto total_volume = 0.0;
        double min_vol {std::numeric_limits<double>::infinity()};
        double max_vol {0.0};
        for (const auto& cell : prism_mesh.cells()) {
            auto cell_vol = cell.volume();
            if (cell_vol < min_vol) {
                min_vol = cell_vol;
            }

            if (cell_vol > max_vol) {
                max_vol = cell_vol;
            }
            total_volume += cell_vol;
        }
        fmt::print("Total volume of mesh: {:.2f} L^3\n", total_volume);
        fmt::print("Max. cell volume = {} L^3\n", max_vol);
        fmt::print("Min. cell volume = {} L^3\n", min_vol);

        // calculate surface area of mesh, using boundary faces area
        auto total_surface_area =
            std::accumulate(prism_mesh.faces().begin(), prism_mesh.faces().end(), 0.0,
                            [](auto acc, const auto& face) { return acc + face.area(); });
        fmt::print("Total surface area of mesh: {:.2f} L^2\n\n", total_surface_area);

        // print number of faces with zero owners
        auto n_zero_owner_faces =
            std::count_if(prism_mesh.faces().begin(), prism_mesh.faces().end(),
                          [](const auto& face) { return !face.has_owner(); });

        // print number of faces with zero owners, and add in green (Ok) or red (Error) at the end of line
        fmt::print("Number of faces with zero owners: {} ", n_zero_owner_faces);
        if (n_zero_owner_faces == 0) {
            fmt::print(fg(fmt::color::green), "(Ok)\n\n");
        } else {
            fmt::print(fg(fmt::color::red), "(Error)\n\n");
        }

        // print boundary conditions
        fmt::print("Mesh boundary patches:\n");
        fmt::print("----------------------\n");
        for (const auto& bc : prism_mesh.boundary_conditions()) {
            fmt::print("  - {} (type: {})\n", bc.name(), boundary_type_to_string(bc.type()));
        }

        // calculate min and max non-orthogonality
        std::vector<double> non_ortho;
        for (const auto& face : prism_mesh.faces()) {
            if (face.neighbor()) {
                non_ortho.push_back(prism_mesh.non_ortho(face));
            }
        }

        auto [min_non_orth, max_non_orth] = min_max(non_ortho);
        fmt::print("\n");
        fmt::print("Max. non-orthogonality = {}\n", max_non_orth);
        fmt::print("Min. non-orthogonality = {}\n", min_non_orth);


    } catch (const std::exception& e) {
        prism::error(e.what());
        return;
    }
}

auto main(int argc, char* argv[]) -> int {
    prism::print_header();

    try {
        cxxopts::Options options("pmesh", "A mesh utility for PMesh files");

        options.set_width(100).set_tab_expansion().add_options()
            // check input mesh and print quality stats
            ("c,check", "Check mesh file", cxxopts::value<std::string>())
            // create new boundary conditions file for the given input mesh, if not found
            ("n,new-boundary", "Create new boundary conditions file for the given input mesh",
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
