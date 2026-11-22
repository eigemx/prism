#pragma once

#include "solver.h"

#include <Eigen/SparseCholesky>

namespace prism::solver {

/** @brief Direct sparse LDL^T solver for the symmetric pressure-correction system.
 *
 * The pressure matrix is symmetric (positive definite after the reference cell is pinned),
 * and for the mesh sizes of interest a sparse Cholesky factorization is fast, exact, and
 * robust. The symbolic analysis only depends on the fixed mesh pattern, so it is reused
 * across factorizations.
 */
class DirectLDLT : public ISolver {
  public:
    void prepare(const SparseMatrix& A) override;
    auto step(const VectorXd& x, const VectorXd& b) -> VectorXd override;

  private:
    Eigen::SimplicialLDLT<SparseMatrix> _ldlt;
    bool _analyzed {false};
    Eigen::Index _n {0};
    Eigen::Index _nnz {0};
};

} // namespace prism::solver
