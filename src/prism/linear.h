#pragma once

#include "types.h"

namespace prism {
// A linear system of equations of the form Ax = b
// where A is a sparse matrix and b is a vector
class LinearSystem {
  public:
    LinearSystem(std::size_t n_cells) : A(n_cells, n_cells), b(n_cells) { b.setZero(); }

    inline auto matrix() const -> const SparseMatrix& { return A; }
    inline auto matrix() -> SparseMatrix& { return A; }
    inline auto matrix(std::size_t i, std::size_t j) const -> double { return A.coeff(i, j); }
    inline auto matrix(std::size_t i, std::size_t j) -> double& { return A.coeffRef(i, j); }

    inline auto rhs() const -> const VectorXd& { return b; }
    inline auto rhs() -> VectorXd& { return b; }
    inline auto rhs(std::size_t i) const -> double { return b[i]; }
    inline auto rhs(std::size_t i) -> double& { return b[i]; }

  private:
    SparseMatrix A;
    VectorXd b;
};

} // namespace prism
