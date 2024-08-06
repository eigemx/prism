#pragma once

#include <prism/types.h>

#include <vector>

namespace prism {
// A linear system of equations of the form Ax = b
class LinearSystem {
  public:
    LinearSystem() = delete;
    LinearSystem(const LinearSystem& other) = default;
    LinearSystem(LinearSystem&& other) noexcept = default;
    auto operator=(const LinearSystem& other) -> LinearSystem& = default;
    auto operator=(LinearSystem&& other) noexcept -> LinearSystem& = default;
    virtual ~LinearSystem() = default;

    LinearSystem(std::size_t n_cells, bool need_matrix = true);

    // setters and getters for matrix A in the right hand side of the linear system
    auto inline matrix() const -> const SparseMatrix& { return _A; }
    auto inline matrix() -> SparseMatrix& { return _A; }

    // insert() updates the matrix coefficient at (i, j), with value v
    // the value v does not overwrite existing value at that location but gets
    // accumulated. inser() does not actually insert directly into matrix A,
    // but actually inserts to a vector of triplets `_triplets` until collect()
    // is called, which does the actual insertion to the matrix.
    void insert(std::size_t i, std::size_t j, double v);

    // collect() should be called after inserting all the matrix coefficients
    void collect();

    // getters and setters for vector `b` in the left hand side of the equation
    auto inline rhs() const -> const VectorXd& { return _b; }
    auto inline rhs() -> VectorXd& { return _b; }

    // getters and setters for the ith-row in the rhs vector.
    auto inline rhs(std::size_t i) const -> double { return _b[i]; }
    auto inline rhs(std::size_t i) -> double& { return _b[i]; }

  private:
    SparseMatrix _A;
    VectorXd _b;
    std::vector<Eigen::Triplet<double>> _triplets;
};

} // namespace prism
