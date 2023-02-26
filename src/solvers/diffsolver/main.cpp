#include <prism/equation.h>
#include <prism/fvscheme.h>
#include <prism/mesh/pmesh.h>
#include <prism/mesh/unv.h>
#include <prism/print.h>
#include <prism/visitor.h>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        prism::error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file = args[1];
    auto mesh = prism::mesh::UnvToPMesh(unv_file).to_pmesh();
    auto eqn = prism::ConservedScalarSteadyEquation("T", mesh);
    auto diffusion_scheme =
        prism::LinearDiffusionScheme(1.0, mesh, eqn.coeff_matrix(), eqn.lhs_vector());

    auto visitor = prism::CellVisitor(mesh, eqn);
    visitor.add_scheme(diffusion_scheme);
    visitor.visit_cells();
    visitor.visit_cells();
}
