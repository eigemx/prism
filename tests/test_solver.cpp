#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>

#include "test_utils.h"

using namespace prism;
using namespace prism::test;

TEST_CASE("All solvers converge to known solution of linear system", "[solver]") {
    constexpr size_t n = 10;
    prism::VectorXd x_analytical(n);
    x_analytical << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    prism::SparseMatrix A;

    // clang-format off
    Eigen::MatrixXd dense(n, n);
        dense <<  6, -1,  0,  0,  0,  0,  0,  0,  0,  0,
                 -1,  6, -1,  0,  0,  0,  0,  0,  0,  0,
                  0, -1,  6, -1,  0,  0,  0,  0,  0,  0,
                  0,  0, -1,  6, -1,  0,  0,  0,  0,  0,
                  0,  0,  0, -1,  6, -1,  0,  0,  0,  0,
                  0,  0,  0,  0, -1,  6, -1,  0,  0,  0,
                  0,  0,  0,  0,  0, -1,  6, -1,  0,  0,
                  0,  0,  0,  0,  0,  0, -1,  6, -1,  0,
                  0,  0,  0,  0,  0,  0,  0, -1,  6, -1,
                  0,  0,  0,  0,  0,  0,  0,  0, -1,  6;
        A = dense.sparseView();
    // clang-format on


    prism::VectorXd b = A * x_analytical;

    SECTION("BiCGSTAB") {
        prism::solver::BiCGSTAB solver;
        prism::VectorXd x = prism::VectorXd::Zero(n);

        solver.prepare(A);
        for (size_t iter = 0; iter < 20; ++iter) {
            x = solver.step(x, b);
        }

        double error = l2NormRel(x, x_analytical);
        REQUIRE(error < 1e-4);
    }

    SECTION("Jacobi") {
        prism::solver::Jacobi solver;
        prism::VectorXd x = prism::VectorXd::Zero(n);

        solver.prepare(A);
        for (size_t iter = 0; iter < 20; ++iter) {
            x = solver.step(x, b);
        }

        double error = l2NormRel(x, x_analytical);
        REQUIRE(error < 1e-4);
    }
}