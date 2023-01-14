#include "mesh/boundary.h"
#include "mesh/from_unv.h"
#include "print.h"


int main() {
    try {
        // load unv mesh file in "./test/diffusion_cube_2D/mesh.unv"
        auto unv_mesh = prism::mesh::UnvToPMesh("./test/slice_2d/slice_2d.unv");
        unv_mesh.report_mesh_stats();

        // convert unv mesh to prism mesh
        auto prism_mesh = unv_mesh.to_pmesh();

        // print boundary conditions
        for (const auto& bc : prism_mesh.boundary_conditions()) {
            fmt::print("Boundary condition: {}\n", bc.name());
        }
    } catch (const std::exception& e) {
        prism::error(e.what());
        return -1;
    }
}
