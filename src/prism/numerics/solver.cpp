#include "solver.h"

namespace prism::solver {

namespace detail {
auto residual(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> double {
    /** @brief Computes the scaled residual of the linear system Ax = b.
     *
     * The scaled residual is defined as:
     * |Ax - b| / max(|A.diagonal() * x|)
     * This is based on equations (14.33) and (14.34) from Moukalled et al. (2015).
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
} // namespace detail

IterationData::IterationData(std::size_t iteration,
                             double initial_residual,
                             double final_residual)
    : _iteration(iteration),
      _initial_residual(initial_residual),
      _final_residual(final_residual) {}

void IterationData::setAsConverged() noexcept {
    _converged = true;
}

auto IterationData::atConvergence() const noexcept -> bool {
    return _converged;
}

auto IterationData::iteration() const noexcept -> std::size_t {
    return _iteration;
}

auto IterationData::initialResidual() const noexcept -> double {
    return _initial_residual;
}

auto IterationData::finalResidual() const noexcept -> double {
    return _final_residual;
}

} // namespace prism::solver
