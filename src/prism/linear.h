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

    RHSProvider(std::size_t n_cells);

    // getters and setters for vector `b` in the left hand side of the equation
    auto inline rhs() const -> const VectorXd& { return _b; }
    auto inline rhs() -> VectorXd& { return _b; }

    // getters and setters for the ith-row in the rhs vector.
    auto inline rhs(std::size_t i) const -> double { return _b[i]; }
    auto inline rhs(std::size_t i) -> double& { return _b[i]; }

  private:
    VectorXd _b;
};

// A linear system of equations of the form Ax = b
class LinearSystem : public RHSProvider {
  public:
    LinearSystem(std::size_t n_cells);

    // setters and getters for matrix A in the right hand side of the linear system
    auto inline matrix() const -> const SparseMatrix& { return _A; }
    auto inline matrix() -> SparseMatrix& { return _A; }

    // insert() updates the matrix coefficient at (i, j), with value v the value v does not
    // overwrite existing value at that location but gets accumulated. insert() does not actually
    // insert directly into matrix A, but inserts to a vector of triplets `_triplets` until
    // collect() is called, which does the actual insertion to the matrix.
    void insert(std::size_t i, std::size_t j, double v);

    // collect() should be called after inserting all the matrix coefficients
    void collect();

  private:
    SparseMatrix _A;
    std::vector<Eigen::Triplet<double>> _triplets;
};

} // namespace prism
