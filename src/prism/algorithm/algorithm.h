#pragma once

#include <span>

#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"

namespace prism::algo {

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
