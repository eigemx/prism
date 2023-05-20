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

        // iterate over the rows of A
        for (std::size_t j = 0; j < n_cells; j++) {
            // dot product the j-th row of A with the solution vector phi using Eigen dot
            auto row_dot_phi = A.row(j).dot(phi.data());

            // subtract the a_jj * phi_j term from the dot product
            auto a_jj = A.coeff(j, j);
            auto& phi_j = phi[j];
            row_dot_phi -= a_jj * phi_j;

            // subtract the dot product from the right hand side vector b
            // and divide by the diagonal element of the j-th row of A
            phi_j = (b[j] - row_dot_phi) / a_jj;
        }


        print("Completed iteration number: {} - residuals norm = {}\n", i, res_norm);

        // zero out the right hand side vector b, so that it can be recalculated
        // in the next iteration using non-orthogonal corrections (if any)
        b.setZero();
    }
}
} // namespace prism::solver