#include "solver.h"

#include <Eigen/IterativeLinearSolvers>

#include "../print.h"
#include "./relax.h"

namespace prism::solver {
void BiCGSTAB::solve(Equation& eqn, std::size_t n_iter, double eps, double lambda) {
    auto n_cells = eqn.scalar_field().data().size();
    auto res = VectorXd(n_cells);

    const auto& A = eqn.coeff_matrix();
    auto& phi = eqn.scalar_field();
    auto& phi_old = eqn.scalar_field_old();
    auto& b = eqn.rhs_vector();

    auto rx = ImplicitUnderRelaxation(lambda);

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.update_coeffs();
        rx.relax(eqn);

        // calculate the norm of the residuals
        res = (A * phi.data()) - b;
        auto res_norm = res.norm();

        // check for convergence
        if (res_norm < eps) {
            fmt::print("Converged after {} iterations\n", i);
            fmt::print("Residual: {}\n", res_norm);
            break;
        }

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
        bicg.compute(A);

        phi_old.data() = phi.data();
        phi.data() -= bicg.solve(res);

        // print the norm of the residuals
        fmt::print("Iteration: {}, Residual: {}\n", i, res_norm);

        // zero out the left & right hand side vector, for the next iteration
        eqn.zero_out_coeffs();
    }
}
} // namespace prism::solver