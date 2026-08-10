#include "jacobi.h"

namespace prism::solver {

void Jacobi::prepare(const SparseMatrix& A) {
    _A = &A;
    _D_inv = A.diagonal().cwiseInverse();
}

auto Jacobi::step(const VectorXd& x, const VectorXd& b) -> VectorXd {
    VectorXd r = b - ((*_A) * x);
    return x + (_D_inv.asDiagonal() * r);
}

} // namespace prism::solver
