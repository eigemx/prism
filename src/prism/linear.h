#pragma once

#include "types.h"

namespace prism {
// A linear system of equations of the form Ax = b
class LinearSystem {
  public:
    LinearSystem(std::size_t n_cells) : A(n_cells, n_cells), b(n_cells) { b.setZero(); }

    auto inline matrix() const -> const SparseMatrix& { return A; }
    auto inline matrix() -> SparseMatrix& { return A; }
    auto inline matrix(std::size_t i, std::size_t j) const -> double { return A.coeff(i, j); }
    auto inline matrix(std::size_t i, std::size_t j) -> double& { return A.coeffRef(i, j); }

    auto inline rhs() const -> const VectorXd& { return b; }
    auto inline rhs() -> VectorXd& { return b; }
    auto inline rhs(std::size_t i) const -> double { return b[i]; }
    auto inline rhs(std::size_t i) -> double& { return b[i]; }

  private:
    SparseMatrix A;
    VectorXd b;
};

} // namespace prism
