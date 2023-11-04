#pragma once

#include <vector>

#include "print.h"
#include "types.h"

namespace prism {
// A linear system of equations of the form Ax = b
class LinearSystem {
  public:
    LinearSystem() = delete;
    LinearSystem(const LinearSystem& other) = default;
    LinearSystem(LinearSystem&& other) = default;
    auto operator=(const LinearSystem& other) -> LinearSystem& = default;
    auto operator=(LinearSystem&& other) -> LinearSystem& = default;
    virtual ~LinearSystem() = default;

    LinearSystem(std::size_t n_cells) : _A(n_cells, n_cells), _b(n_cells) {
        // set the right hand side to zero
        _b.setZero();

        // assume that the mesh is purely tetrahedral, so each cell has 3 neighbors
        // then the number of triplets (i, j, v) is 4 * n_cells
        // this is the minimum number of triplets, but it is not the exact number
        // In the general case of a polyhedral mesh, this number is not correct
        // but can be treated as a warm-up for allocations.
        // TODO: check if this is making any difference, and remove it if not.
        _triplets.reserve(4 * n_cells);
    }

    // setters and getters for matrix A in the right hand side of the linear system
    auto inline matrix() const -> const SparseMatrix& { return _A; }
    auto inline matrix() -> SparseMatrix& { return _A; }

    // insert() updates the matrix coefficient at (i, j), with value v
    // the value v does not overwrite existing value at that location but gets
    // accumulated. inser() does not actually insert directly into matrix A,
    // but actually inserts to a vector of triplets `_triplets` until collect()
    // is called, which does the actual insertion to the matrix.
    auto inline insert(std::size_t i, std::size_t j, double v) -> void {
        _triplets.emplace_back(i, j, v);
    }

    // collect() should be called after inserting all the matrix coefficients
    void inline collect() {
        if (_triplets.empty()) {
            throw std::runtime_error(
                "LinearSystem::collect() was called on an empty "
                "triplet list. This should not happen.");
        }
        _A.setFromTriplets(_triplets.begin(), _triplets.end());
        _triplets.clear();
    }

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
