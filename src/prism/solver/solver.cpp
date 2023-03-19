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

        // start the i-th iteration
        for (std::size_t j = 0; j < n_cells; j++) {
            double sum = 0.0;

            for (SparseMatrix::InnerIterator it(A, j); it; ++it) {
                if (it.row() != it.col()) {
                    sum += it.value() * phi[it.row()];
                }
            }
            phi[j] = (b[j] - sum) / A.coeff(j, j);

            // update the jth residual
            res[j] = std::abs(b[j] - (A.row(j) * phi.data()).sum());
        }

        // check root mean square of residual
        double rms = res.norm() / std::sqrt(n_cells);

        print("Iteration {}: RMS = {}\n", i, rms);

        if (rms < eps) {
            print("Converged in {} iterations\n", i);
            return;
        }
    }
}
} // namespace prism::solver