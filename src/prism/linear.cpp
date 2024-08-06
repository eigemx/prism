#include "linear.h"

namespace prism {
LinearSystem::LinearSystem(std::size_t n_cells, bool need_matrix)
    : _A(n_cells, n_cells), _b(n_cells) {
    // set the right hand side to zero
    _b.setZero();

    if (need_matrix) {
        // assume that the mesh is purely tetrahedral, so each cell has 3 neighbors
        // then the number of triplets (i, j, v) is 4 * n_cells
        // this is the minimum number of triplets, but it is not the exact number
        // In the general case of a polyhedral mesh, this number is not correct
        // but can be treated as a warm-up for allocations.
        // TODO: check if this is making any difference, and remove it if not.
        _triplets.reserve(4 * n_cells);
    }
}

void LinearSystem::insert(std::size_t i, std::size_t j, double v) {
    _triplets.emplace_back(i, j, v);
}

void LinearSystem::collect() {
    if (_triplets.empty()) {
        throw std::runtime_error(
            "LinearSystem::collect() was called on an empty "
            "triplet list. This should not happen.");
    }
    _A.setFromTriplets(_triplets.begin(), _triplets.end());
    _triplets.clear();
}
} // namespace prism