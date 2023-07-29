#pragma once

#include <Eigen/IterativeLinearSolvers>

#include "../equation.h"
#include "../field.h"
#include "../print.h"
#include "../types.h"
#include "relax.h"

namespace prism::solver {
// TODO: make this a template class, with the template parameter being the solver type
// TODO: solve() should return each time an iteraion is finished (with the current residual),
// until convergence. And any consecutive call to solve() after convergence shall have no effect
class SolverBase {
  public:
    SolverBase() = default;
    SolverBase(const SolverBase& s) = delete;
    SolverBase(SolverBase&& s) = delete;
    auto operator=(SolverBase&& s) -> SolverBase& = delete;
    auto operator=(const SolverBase& s) -> SolverBase& = delete;
    virtual ~SolverBase() = default;

    virtual void solve(Equation& eq, std::size_t n_iter, double eps, double lambda) = 0;
};

template <typename Relaxer = ImplicitUnderRelaxation>
class BiCGSTAB : public SolverBase {
  public:
    void solve(Equation& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-4,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) override;
};

template <typename Relaxer>
void BiCGSTAB<Relaxer>::solve(Equation& eqn, std::size_t n_iter, double eps, double lambda) {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();

    auto& phi = eqn.field();
    auto& phi_prev = eqn.field_prev_iter();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    auto rx = Relaxer(lambda);

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.update_coeffs();
        rx.pre_relax(eqn);

        // calculate the residuals and its norm
        auto res = (A * phi.data()) - b;
        auto res_norm = res.norm();

        // check for convergence
        if (res_norm < eps) {
            fmt::println("Converged after {} iterations", i);
            fmt::println("Residual: {}", res_norm);
            break;
        }

        bicg.compute(A);

        phi_prev.data() = phi.data();
        phi.data() -= bicg.solve(res);

        rx.post_relax(eqn);

        fmt::println("Iteration: {}, Residual: {}", i, res_norm);

        // zero out the left & right hand side vector, for the next iteration
        eqn.zero_out_coeffs();
    }
}

} // namespace prism::solver