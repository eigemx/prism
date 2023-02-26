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
    using namespace prism;

    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        prism::error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file = args[1];

    // read mesh
    auto mesh = mesh::UnvToPMesh(unv_file).to_pmesh();

    // setup steady diffusion equation for temperature
    auto eqn = ConservedScalarSteadyEquation("T", mesh);

    // linear diffusion interpolation scheme
    auto diff_scheme = LinearDiffusionScheme(1.0, mesh, eqn.coeff_matrix(), eqn.lhs_vector());

    // prepare cell visitor and connect to the scheme
    auto visitor = CellVisitor(mesh, eqn);
    visitor.add_scheme(diff_scheme);

    // visit cells and adjust coefficients matrix A and rhs vector b -> A * T = b
    visitor.visit_cells();
}
