#include <prism/core.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    print("diffsolver - A steady state temperature diffusion solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    print("Loading mesh file {}...", unv_file_name);
    auto mesh = mesh::UnvToPMesh(unv_file_name).to_pmesh();
    print("Okay.\n");

    auto T = ScalarField("temperature", mesh, 300.0);

    auto diff = diffusion::Linear(1.0, T);

    auto eqn = Equation(T, {&diff});

    auto solver = solver::GaussSeidel();
    solver.solve(eqn, 10000, 1e-6);

    prism::export_field(eqn.scalar_field(), "sol.vtu");

    /*for (std::size_t i = 0; i < 3; ++i) {
        eqn.update_coeffs();

        print("Completed iteration number: {}\n", i);
    }*/

    //const auto& result = eqn.scalar_field().data();
    //std::ofstream file("temperature.csv");
    //file << result;
    //file.close();


    // export right hand side vector b to csv file `vector.csv`
    /*auto b = eqn.rhs_vector();
    std::ofstream vector_file("vector.csv");
    vector_file << b;
    vector_file.close();

    // export coefficients matrix A to csv file `matrix.csv`
    SparseMatrix& A = eqn.coeff_matrix();
    std::ofstream matrix_file("matrix.csv");
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
        info("Matrix A is solvable.");
    }*/
}
