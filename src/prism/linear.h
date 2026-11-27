#pragma once

#include <prism/types.h>

#include <vector>

namespace prism {

class RHSProvider {
  public:
    RHSProvider() = delete;
    RHSProvider(const RHSProvider& other) = default;
    RHSProvider(RHSProvider&& other) noexcept = default;
    auto operator=(const RHSProvider& other) -> RHSProvider& = default;
    auto operator=(RHSProvider&& other) noexcept -> RHSProvider& = default;
    virtual ~RHSProvider() = default;

    RHSProvider(size_t n_cells);

    // getters and setters for vector `b` in the left hand side of the equation
    auto rhs() const -> const VectorXd& { return _b; }
    auto rhs() -> VectorXd& { return _b; }

    // getters and setters for the ith-row in the rhs vector.
    auto rhs(size_t i) const -> f64 { return _b[static_cast<Eigen::Index>(i)]; }
    auto rhs(size_t i) -> f64& { return _b[static_cast<Eigen::Index>(i)]; }

  private:
    VectorXd _b;
};

// A linear system of equations of the form Ax = b
class LinearSystem : public RHSProvider {
  public:
    LinearSystem(size_t n_cells);

    // setters and getters for matrix A in the right hand side of the linear system
    auto matrix() const -> const SparseMatrix& { return _A; }
    auto matrix() -> SparseMatrix& { return _A; }

    // Adopts a precomputed sparsity pattern (see mesh::PMesh::equationPattern()).
    // Insertions are accumulated directly into the matrix values, so the pattern
    // is never rebuilt and the assembly order does not matter.
    void setPattern(const SparseMatrix& pattern);

    // insert() accumulates the coefficient at (i, j). With a pattern set, (i, j) must
    // exist in it; otherwise coefficients are buffered as triplets until collect().
    void insert(size_t i, size_t j, f64 v);

    // collect() should be called after inserting all the matrix coefficients.
    void collect();

    // Clears the coefficient values and RHS while keeping the sparsity pattern.
    void resetValues();

    // Accumulates the values of another system that shares this system's pattern.
    void accumulateFrom(const LinearSystem& other);

  private:
    SparseMatrix _A;
    std::vector<Eigen::Triplet<f64>> _triplets;

    // Value slot of every matrix entry, grouped by row and sorted by column.
    std::vector<std::vector<Pair<size_t, Eigen::Index>>> _row_slots;
    bool _has_pattern {false};
};

} // namespace prism
