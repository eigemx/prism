#include <filesystem>

#include "mesh/boundary.h"
#include "mesh/from_unv.h"
#include "print.h"


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

auto main() -> int {
    try {
        // load unv mesh file in "./test/diffusion_cube_2D/mesh.unv"
        auto mesh_file = std::filesystem::path {"./test/slice_2d/mesh.unv"};

        fmt::print("Loading mesh file: ");
        fmt::print(fg(fmt::color::dark_cyan), "{}\n", mesh_file.string());

        auto unv_mesh = prism::mesh::UnvToPMesh(mesh_file);
        //unv_mesh.report_mesh_stats();

        // convert unv mesh to prism mesh
        auto prism_mesh = unv_mesh.to_pmesh();

        // print number of cells and faces
        fmt::print("Number of cells: {}\n", prism_mesh.cells().size());
        fmt::print("Number of faces: {}\n\n", prism_mesh.faces().size());

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
        for (const auto& bc : prism_mesh.boundary_conditions()) {
            fmt::print("  - {} (type: {})\n", bc.name(), boundary_type_to_string(bc.type()));
        }
    } catch (const std::exception& e) {
        prism::error(e.what());
        return -1;
    }
}
