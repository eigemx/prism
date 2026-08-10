#pragma once

#include "solver.h"

namespace prism::solver {

class BiCGSTAB : public ISolver {
  public:
    void prepare(const SparseMatrix& A) override;
    auto step(const VectorXd& x, const VectorXd& b) -> VectorXd override;

  private:
    // Eigen's default IncompleteLUT (fill factor 10) is kept because a zero-fill DiLU
    // preconditioner was too weak: restarted BiCGSTAB stalled on the pressure equation and
    // the cavity PISO case diverged. Reducing the fill factor below 10 also destabilised
    // PISO, so the default is retained as the robust choice.
    Eigen::BiCGSTAB<Eigen::SparseMatrix<f64>, Eigen::IncompleteLUT<f64>> _bicg;
};

} // namespace prism::solver
