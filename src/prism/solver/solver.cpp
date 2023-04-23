#include "solver.h"

#include "../print.h"


namespace prism::solver {
void GaussSeidel::solve(Equation& eqn, std::size_t n_iter, double eps) {
    auto n_cells = eqn.scalar_field().data().size();
    auto res = VectorXd(n_cells);

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.update_coeffs();

        auto& A = eqn.coeff_matrix();
        auto& b = eqn.rhs_vector();
        auto& phi = eqn.scalar_field();

        // solve linear system A * phi = b, and check for convergence
        res = (A * phi.data()) - b;
        auto res_norm = res.norm();

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
    }
}
} // namespace prism::solver