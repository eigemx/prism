#include "solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <iostream>

#include "../print.h"

namespace prism::solver {
void BiCGSTAB::solve(Equation& eqn, std::size_t n_iter, double eps) {
    auto n_cells = eqn.scalar_field().data().size();
    auto res = VectorXd(n_cells);

    const auto& A = eqn.coeff_matrix();
    auto& phi = eqn.scalar_field();
    auto& b = eqn.rhs_vector();

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.update_coeffs();

        // calculate the norm of the residuals
        res = (A * phi.data()) - b;
        auto res_norm = res.norm();

        // check for convergence
        if (res_norm < eps) {
            print("Converged after {} iterations\n", i);
            print("Residual: {}\n", res_norm);
            break;
        }

        //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
        //cg.compute(A);
        //phi.data() -= cg.solve(res);

        // solve using BiCGSTAB
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
        bicg.compute(A);
        phi.data() -= bicg.solve(res);

        // print the norm of the residuals
        print("Iteration: {}, Residual: {}\n", i, res_norm);

        // zero out the left & right hand side vector, for the next iteration
        eqn.zero_out_coeffs();
    }
}
} // namespace prism::solver