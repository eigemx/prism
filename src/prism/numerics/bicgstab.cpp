#include "bicgstab.h"

namespace prism::solver {

void BiCGSTAB::prepare(const SparseMatrix& A) {
    // The outer ISolver driver checks the Moukalled-scaled residual. This inner solve is
    // configured separately; the defaults are an empirical budget for the diagonal
    // preconditioner (single-iteration restarts were stable only with ILUT).
    _bicg.setMaxIterations(static_cast<Eigen::Index>(_options.max_iterations));
    _bicg.setTolerance(_options.relative_tolerance);
    _bicg.compute(A);
}

auto BiCGSTAB::step(const VectorXd& x, const VectorXd& b) -> VectorXd {
    return _bicg.solveWithGuess(b, x);
}

} // namespace prism::solver
