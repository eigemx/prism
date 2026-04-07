#include "solver.h"

#include "prism/log.h"

namespace prism::solver {

namespace {
auto residual(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> double {
    /** @brief Computes the scaled residual of the linear system Ax = b.
     *
     * The scaled residual is defined as:
     * |Ax - b| / max(|A.diagonal() * x|)
     * This is based on equations (14.33) and (14.34) from Moukalled et al. (2016).
     *
     * @param A The sparse matrix representing the linear system.
     * @param x The solution vector.
     * @param b The right-hand side vector.
     * @return The scaled residual value.
     */
    auto ac_phic = A.diagonal().cwiseProduct(x);
    auto res_scaled = ((A * x) - b).cwiseAbs() / (ac_phic.cwiseAbs().maxCoeff() + EPSILON);

    return res_scaled.maxCoeff();
}
} // namespace

auto ISolver::solve(eqn::Transport& eqn, std::size_t n_iter, double eps) -> SolverResult {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();
    auto phi = eqn.field();

    f64 init_res = 0.0;
    f64 current_res = 0.0;
    bool has_converged = false;
    size_t iter = 0;

    eqn.updateCoeffs();
    eqn.relax();

    for (; iter < n_iter; iter++) {
        if (iter == 0) {
            init_res = residual(A, phi->values(), b);
        }

        phi->values() = step(A, phi->values(), b);
        current_res = residual(A, phi->values(), b);

        if (current_res < eps) {
            has_converged = true;
            iter++;
            break;
        }
    }
    eqn.zeroOutCoeffs();

    /// TODO: we shouldn't log::info inside ISolver::solve(), executables should handles this.
    /// Remove.
    log::info("Residuals: Initial = {:.4e} | Final: {:.4e} (nIterations = {})",
              init_res,
              current_res,
              iter);

    SolverResult data(iter, init_res, current_res);
    if (has_converged) {
        data.setAsConverged();
    }
    return data;
}

SolverResult::SolverResult(std::size_t iteration, double initial_residual, double final_residual)
    : _iteration(iteration), _initial_residual(initial_residual), _final_residual(final_residual) {}

void SolverResult::setAsConverged() noexcept {
    _converged = true;
}

auto SolverResult::hasConverged() const noexcept -> bool {
    return _converged;
}

auto SolverResult::iteration() const noexcept -> std::size_t {
    return _iteration;
}

auto SolverResult::initialResidual() const noexcept -> double {
    return _initial_residual;
}

auto SolverResult::finalResidual() const noexcept -> double {
    return _final_residual;
}

} // namespace prism::solver