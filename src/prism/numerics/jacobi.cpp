#include "jacobi.h"

namespace prism::solver {

auto Jacobi::step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd {
    VectorXd D_inv = A.diagonal().cwiseInverse();
    VectorXd r = b - (A * x);
    return x + (D_inv.asDiagonal() * r);
}

} // namespace prism::solver