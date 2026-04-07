#pragma once

#include "solver.h"

namespace prism::solver {

class BiCGSTAB : public ISolver {
  public:
    auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd override;

  private:
    Eigen::BiCGSTAB<Eigen::SparseMatrix<f64>, Eigen::IncompleteLUT<f64>> _bicg;
};

} // namespace prism::solver