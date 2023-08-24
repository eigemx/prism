#pragma once

#include <vector>

#include "types.h"

namespace prism {
// A linear system of equations of the form Ax = b
class LinearSystem {
  public:
    LinearSystem(std::size_t n_cells) : _A(n_cells, n_cells), _b(n_cells) {
        _b.setZero();

        // assume that the mesh is purely tetrahedral, so each cell has 3 neighbors
        // then the number of triplets (i, j, v) is 4 * n_cells
        // this is the minimum number of triplets, but it is not the exact number
        _triplets.reserve(4 * n_cells);
    }

    auto inline matrix() const -> const SparseMatrix& { return _A; }
    auto inline matrix() -> SparseMatrix& { return _A; }

    auto inline insert(std::size_t i, std::size_t j, double v) -> void {
        _triplets.emplace_back(i, j, v);
    }

    void inline collect() {
        _A.setFromTriplets(_triplets.begin(), _triplets.end());
        _triplets.clear();
    }

    auto inline rhs() const -> const VectorXd& { return _b; }
    auto inline rhs() -> VectorXd& { return _b; }
    auto inline rhs(std::size_t i) const -> double { return _b[i]; }
    auto inline rhs(std::size_t i) -> double& { return _b[i]; }

  private:
    SparseMatrix _A;
    VectorXd _b;
    std::vector<Eigen::Triplet<double>> _triplets;
};

} // namespace prism
