#pragma once

#include <Eigen/IterativeLinearSolvers>

#include "solver.h"

namespace prism::solver {

/** @brief Inner BiCGSTAB configuration.
 *
 * These values configure Eigen's Krylov solve inside step(). They are separate from
 * SolverOptions, which controls the outer ISolver driver and its scaled-residual check.
 */
struct BiCGSTABOptions {
    // Maximum inner Krylov iterations per step().
    size_t max_iterations {100};

    // Inner stopping tolerance: Eigen's relative residual |Ax - b| / |b|.
    f64 relative_tolerance {1e-8};
};

class BiCGSTAB : public ISolver {
  public:
    explicit BiCGSTAB(BiCGSTABOptions options = {}) : _options(options) {}

    void prepare(const SparseMatrix& A) override;
    auto step(const VectorXd& x, const VectorXd& b) -> VectorXd override;

  private:
    BiCGSTABOptions _options;

    // Momentum is nonsymmetric but diagonally dominant after under-relaxation. A diagonal
    // preconditioner has O(n) setup and a cheap apply.
    Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<f64>> _bicg;
};

} // namespace prism::solver
