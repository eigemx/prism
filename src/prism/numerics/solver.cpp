#include "solver.h"

#include "prism/log.h"

namespace prism::solver {

namespace {
auto residual(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> f64 {
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

auto ISolver::solve(eqn::Transport& eqn, const SolverOptions& options) -> SolverResult {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();
    auto phi = eqn.field();

    f64 init_res = 0.0;
    f64 current_res = 0.0;
    bool has_converged = false;
    size_t iter = 0;

    eqn.updateCoeffs();
    eqn.relax();

    // OpenFOAM's setReference: pin the solution at one cell to break the null
    // space of all-Neumann (closed-domain) pressure systems. Doubling the
    // diagonal removes the constant vector from the null space since the row
    // sum is no longer zero. The RHS adjustment b[ref] += diag * ref_value
    // anchors the solution to the desired reference value.
    if (options.reference_cell.has_value() &&
        *options.reference_cell < size_t(eqn.matrix().rows())) {
        auto& mat = eqn.matrix();
        auto& vec = eqn.rhs();
        auto idx = static_cast<Eigen::Index>(*options.reference_cell);
        auto diag = mat.coeff(idx, idx);
        mat.coeffRef(idx, idx) += diag;
        vec[idx] += diag * options.reference_value;
    }

    // Factorize the (relaxed, reference-pinned) matrix once; step() then reuses the
    // factorization for every iteration instead of recomputing it each time. Convergence is
    // driven by the scaled residual check below rather than by any internal solver tolerance.
    if (options.reuse_prepared) {
        reattach(A);
    } else {
        prepare(A);
    }

    for (; iter < options.max_iterations; iter++) {
        if (iter == 0) {
            init_res = residual(A, phi->values(), b);
        }

        phi->values() = step(phi->values(), b);
        current_res = residual(A, phi->values(), b);

        if (current_res < options.tolerance) {
            has_converged = true;
            iter++;
            break;
        }
    }
    if (options.clear_coeffs) {
        eqn.zeroOutCoeffs();
    }

    SolverResult data(iter, init_res, current_res);
    if (has_converged) {
        data.setAsConverged();
    }
    return data;
}

SolverResult::SolverResult(size_t iteration, f64 initial_residual, f64 final_residual)
    : _iteration(iteration), _initial_residual(initial_residual), _final_residual(final_residual) {}

void SolverResult::setAsConverged() noexcept {
    _converged = true;
}

auto SolverResult::hasConverged() const noexcept -> bool {
    return _converged;
}

auto SolverResult::iteration() const noexcept -> size_t {
    return _iteration;
}

auto SolverResult::initialResidual() const noexcept -> f64 {
    return _initial_residual;
}

auto SolverResult::finalResidual() const noexcept -> f64 {
    return _final_residual;
}

} // namespace prism::solver