#pragma once

#include <span>

#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"

namespace prism::algo {

struct StepResult {
    auto hasConverged() const -> bool;
    auto hasDiverged() const -> bool;
    auto momentumResiduals() const -> std::array<double, 3>;
    auto pressureResidual() const -> double;

  private:
    bool _has_converged = false;
    bool _has_diverged = false;
    std::array<double, 3> _momentum_residuals = {0., 0., 0.};
    double _pressure_residual = 0.;
};

class IPressureLinked {
  public:
    IPressureLinked() = default;
    IPressureLinked(const IPressureLinked&) = default;
    IPressureLinked(IPressureLinked&&) = delete;
    auto operator=(const IPressureLinked&) -> IPressureLinked& = default;
    auto operator=(IPressureLinked&&) -> IPressureLinked& = delete;
    virtual ~IPressureLinked() = default;

    virtual void step(std::span<eqn::Momentum*> momentum_predictors,
                      field::Velocity& U,
                      field::Velocity& mdot,
                      field::Pressure& p) = 0;
};
} // namespace prism::algo
