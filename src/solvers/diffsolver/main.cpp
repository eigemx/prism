#include <prism/core.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    /*using namespace prism;

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
    std::size_t n_row_non_zero = 0;

    for (int i = 0; i < A.rows(); i++) {
        if ((b[i] > 1e-8) || (b[i] < -1e-8)) {
            continue;
        }
        double row_sum = 0.0;
        for (int j = 0; j < A.cols(); j++) {
            row_sum += A.coeff(i, j);
        }

        if (std::abs(row_sum) > 1e-8) {
            n_row_non_zero++;
        }
    }

    if (n_row_non_zero > 0) {
        error(format("Matrix A has {} rows that do not sum to zero", n_row_non_zero));
    } else {
        info("Matrix A is solvable\n");
    }*/

    using namespace prism;

    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto mesh = mesh::UnvToPMesh(unv_file_name).to_pmesh();

    auto T = ScalarField("T", mesh, 300.0);
    auto V = VectorField("V", mesh, Vector3d {1.0, 2.0, 3.0});

    auto V_x = V.x();

    // update values of V_x by adding 10.0 to each elelment, using Eigen::VectorXd methods
    for (auto& val : V_x.data()) {
        val += 10.0;
    }

    V_x.update_parent_vec_field();

    std::cout << V.data() << std::endl;
}
