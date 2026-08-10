#pragma once

#include "solver.h"

namespace prism::solver {

class Jacobi : public ISolver {
  public:
    void prepare(const SparseMatrix& A) override;
    auto step(const VectorXd& x, const VectorXd& b) -> VectorXd override;

  private:
    const SparseMatrix* _A {nullptr};
    VectorXd _D_inv;
};

} // namespace prism::solver
