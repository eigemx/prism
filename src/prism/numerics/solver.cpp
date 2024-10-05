#include "solver.h"

namespace prism::solver {

IterationData::IterationData(std::size_t iteration,
                             double initial_residual,
                             double final_residual)
    : _iteration(iteration),
      _initial_residual(initial_residual),
      _final_residual(final_residual) {}

void IterationData::setAsConverged() noexcept {
    _converged = true;
}

auto IterationData::atConvergence() const noexcept -> bool {
    return _converged;
}

auto IterationData::iteration() const noexcept -> std::size_t {
    return _iteration;
}

auto IterationData::initialResidual() const noexcept -> double {
    return _initial_residual;
}

auto IterationData::finalResidual() const noexcept -> double {
    return _final_residual;
}

} // namespace prism::solver
