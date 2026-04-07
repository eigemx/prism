#include "bicgstab.h"

#include <Eigen/IterativeLinearSolvers>

namespace prism::solver {

auto BiCGSTAB::step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd {
    return _bicg.compute(A).solveWithGuess(b, x);
}

} // namespace prism::solver