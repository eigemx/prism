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
    auto eqn = SteadyConservedScalar("T", mesh);

    // Add diffusion term (∇.k∇T) (linear diffusion interpolation scheme)
    auto diff_scheme = diffusion::Linear(1.0);

    // add explicit gradient scheme (to calculate ∇T in non-orthogonal corrections)
    auto grad_scheme = gradient::GreenGauss();

    // prepare cell visitor and connect to the scheme
    //eqn.add_scheme(diff_scheme);
    eqn.add_scheme(diffusion::Linear(1.0));

    eqn.update_coeffs();

    // export right hand side vector b to csv file `vector.csv`
    auto b = eqn.lhs_vector();
    std::ofstream vector_file("vector_u.csv");
    vector_file << b;
    vector_file.close();

    // export coefficients matrix A to csv file `matrix.csv`
    SparseMatrix& A = eqn.coeff_matrix();
    std::ofstream matrix_file("matrix_u.csv");
    matrix_file << A;
    matrix_file.close();

    // check that for each row of matrix A, the row sums to zero
    // (this is a requirement for the linear system to be solvable)
    for (int i = 0; i < A.rows(); i++) {
        double row_sum = 0.0;
        for (int j = 0; j < A.cols(); j++) {
            row_sum += A.coeff(i, j);
        }

        row_sum -= b[i];
        if (std::abs(row_sum) > 1e-8) {
            error(format("Row {} of matrix A does not sum to zero", i));
        }
    }
}
