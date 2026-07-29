#pragma once

#include <span>

#include "algorithm.h"
#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"
#include "prism/field/velocity.h"
#include "prism/report.h"
#include "prism/types.h"

namespace prism::algo {
struct PRIMEControls {
    // Number of non-orthogonal correctors for pressure equation
    size_t non_ortho_correctors = 2;

    // Under-relaxation factor for pressure field
    f64 pressure_urf = 0.3;

    // Pressure equation inner solve() loop max. iterations count
    size_t pressure_max_iter = 5;

    // Minimum residual for pressure equation
    f64 pressure_residual = 1e-7;

    // Reference cell for singular (all-Neumann) pressure systems.
    Optional<size_t> p_ref_cell = NullOption;
    f64 p_ref_value = 0.0;
};


class PRIME : public IPressureLinked {
  public:
    PRIME(PRIMEControls controls);
    auto step(std::span<eqn::Momentum*> momentum_predictors,
              SharedPtr<field::Velocity>& U,
              SharedPtr<field::Velocity>& mdot,
              SharedPtr<field::Pressure>& p) -> std::vector<report::Entry> override;

  private:
    PRIMEControls _controls;
};

void solveExplicitMomentum(std::span<eqn::Momentum*> momentum_predictors);

} // namespace prism::algo
