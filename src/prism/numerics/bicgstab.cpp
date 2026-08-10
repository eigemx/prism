#include "bicgstab.h"

#include <Eigen/IterativeLinearSolvers>

namespace prism::solver {

void BiCGSTAB::prepare(const SparseMatrix& A) {
    // step() performs exactly one BiCGSTAB iteration; the ISolver driver checks its own
    // (Moukalled-scaled) residual between iterations and stops once it is below the
    // configured tolerance. Restarting after every single iteration avoids the BiCGSTAB
    // breakdown (rho -> 0) that longer uninterrupted runs hit on advection-dominated
    // systems, while the driver's residual remains the convergence criterion (so the
    // solver does not over-solve to Eigen's machine-precision default).
    _bicg.setMaxIterations(1);
    _bicg.setTolerance(0.0);
    _bicg.compute(A);
}

auto BiCGSTAB::step(const VectorXd& x, const VectorXd& b) -> VectorXd {
    return _bicg.solveWithGuess(b, x);
}

} // namespace prism::solver
