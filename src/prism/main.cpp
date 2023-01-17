#include "mesh/boundary.h"
#include "mesh/from_unv.h"
#include "print.h"


auto main() -> int {
    try {
        // load unv mesh file in "./test/diffusion_cube_2D/mesh.unv"
        auto unv_mesh = prism::mesh::UnvToPMesh("./test/cylinder/mesh.unv");
        //unv_mesh.report_mesh_stats();

        // convert unv mesh to prism mesh
        auto prism_mesh = unv_mesh.to_pmesh();

        // print number of cells and faces
        fmt::print("Number of cells: {}\n", prism_mesh.cells().size());
        fmt::print("Number of faces: {}\n\n", prism_mesh.faces().size());

        // print boundary conditions
        for (const auto& bc : prism_mesh.boundary_conditions()) {
            fmt::print("Boundary condition: {}\n", bc.name());
        }
    } catch (const std::exception& e) {
        prism::error(e.what());
        return -1;
    }
}
