#include <prism/core.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    print("diffsolver - A steady state temperature diffusion solver\n");

    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto mesh = mesh::UnvToPMesh(unv_file_name).to_pmesh();

    // setup steady diffusion equation for temperature
    // this will create a linear system of equations A * T = b
    auto eqn = ConservedScalarSteadyEquation("T", mesh);

    // Add diffusion term (∇.k∇T) (linear diffusion interpolation scheme)
    auto diff_scheme = diffusion::Linear(1.0, mesh, eqn.coeff_matrix(), eqn.lhs_vector());

    // add explicit gradient scheme (to calculate ∇T in non-orthogonal corrections)
    auto grad_scheme = gradient::GreenGauss(mesh, eqn.coeff_matrix(), eqn.lhs_vector());

    // prepare cell visitor and connect to the scheme
    auto visitor = CellVisitor(mesh, eqn);
    visitor.add_scheme(diff_scheme);
    //visitor.add_scheme(grad_scheme);

    // visit cells and adjust coefficients matrix A and rhs vector b -> A * T = b
    visitor.visit_cells();

    // export coefficients matrix A to csv file `matrix.csv`
    SparseMatrix& A = eqn.coeff_matrix();
    std::ofstream matrix_file("matrix.csv");
    matrix_file << A;
    matrix_file.close();
}
