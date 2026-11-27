#include "linear.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace prism {
RHSProvider::RHSProvider(size_t n_cells) : _b(n_cells) {
    _b.setZero();
}

LinearSystem::LinearSystem(size_t n_cells)
    : RHSProvider(n_cells),
      _A(static_cast<Eigen::Index>(n_cells), static_cast<Eigen::Index>(n_cells)) {}

void LinearSystem::setPattern(const SparseMatrix& pattern) {
    assert(!_has_pattern);
    assert(pattern.rows() == _A.rows() && pattern.cols() == _A.cols());

    _A = pattern;
    _A.coeffs().setZero();
    _triplets.clear();

    // Record the value slot of every entry, grouped by row and ordered by column so
    // insert() can locate a coefficient with lower_bound(). The iterator exposes the
    // real row/column in both storage orders, and appends each row in ascending column
    // order, so no sort is needed.
    _row_slots.assign(_A.rows(), {});
    for (Eigen::Index outer = 0; outer < _A.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator it(_A, outer); it; ++it) {
            const auto slot = static_cast<Eigen::Index>(&it.value() - _A.valuePtr());
            _row_slots[it.row()].emplace_back(static_cast<size_t>(it.col()), slot);
        }
    }

    _has_pattern = true;
}

void LinearSystem::resetValues() {
    if (_has_pattern) {
        _A.coeffs().setZero();
    } else {
        _A.setZero();
        _triplets.clear();
    }
    rhs().setZero();
}

void LinearSystem::accumulateFrom(const LinearSystem& other) {
    assert(_A.rows() == other._A.rows());
    assert(_A.nonZeros() == other._A.nonZeros());
    assert(std::equal(
        _A.outerIndexPtr(), _A.outerIndexPtr() + _A.outerSize() + 1, other._A.outerIndexPtr()));
    assert(std::equal(
        _A.innerIndexPtr(), _A.innerIndexPtr() + _A.nonZeros(), other._A.innerIndexPtr()));

    for (Eigen::Index i = 0; i < _A.nonZeros(); ++i) {
        _A.valuePtr()[i] += other._A.valuePtr()[i];
    }
}

void LinearSystem::insert(size_t i, size_t j, f64 v) {
    if (!_has_pattern) {
        _triplets.emplace_back(i, j, v);
        return;
    }

    const auto& row = _row_slots[i];
    const auto entry = std::lower_bound(
        row.begin(), row.end(), j, [](const auto& lhs, size_t col) { return lhs.first < col; });
    assert(entry != row.end() && entry->first == j);
    _A.valuePtr()[entry->second] += v;
}

void LinearSystem::collect() {
    if (_has_pattern) {
        return;
    }

    if (_triplets.empty()) {
        throw std::runtime_error("LinearSystem::collect() was called on an empty triplet list.");
    }

    _A.setFromTriplets(_triplets.begin(), _triplets.end());
    _triplets.clear();
}
} // namespace prism
