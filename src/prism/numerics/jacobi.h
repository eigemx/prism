#pragma once

#include "solver.h"

namespace prism::solver {

class Jacobi : public ISolver {
  public:
    auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd override;
};

} // namespace prism::solver