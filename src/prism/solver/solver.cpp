#include "solver.h"

#include "../print.h"


namespace prism::solver {
void GaussSeidel::solve(Equation& eqn, std::size_t n_iter, double eps) {
    auto n_cells = eqn.scalar_field().data().size();
    auto res = VectorXd(n_cells);

    const auto& A = eqn.coeff_matrix();
    auto& phi = eqn.scalar_field();
    auto& b = eqn.rhs_vector();

    for (std::size_t i = 0; i < n_iter; i++) {
        // Warning: This corrects non-orthogonality in each iteration.
        eqn.update_coeffs();

        // calculate the norm of the residuals
        res = (A * phi.data()) - b;
        auto res_norm = res.norm();

        // check for convergence
        if (res_norm < eps) {
            print("Converged after {} iterations\n", i);
            break;
        }

        // solve linear system A * phi = b
        for (int j = 0; j < A.rows(); j++) {
            double row_sum = 0.0;
            for (int k = 0; k < A.cols(); k++) {
                if (j != k) {
                    row_sum += A.coeff(j, k) * phi.data()[k];
                }
            }
            phi[j] = (b[j] - row_sum) / A.coeff(j, j);
        }

        print("Completed iteration number: {} - residuals norm = {}\n", i, res_norm);

        // zero out the right hand side vector b, so that it can be recalculated
        // in the next iteration using non-orthogonal corrections (if any)
        b.setZero();
    }
}
} // namespace prism::solver