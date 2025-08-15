#pragma once

#include <cctype>
#include <span>

#include "algorithm.h"
#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"
#include "prism/field/velocity.h"

namespace prism::algo {
struct PRIMEParameters {
    // Number of non-orthogonal correctors for pressure equation
    std::size_t non_ortho_correctors = 2;

    // Under-relaxation factor for pressure field
    double pressure_urf = 0.3;

    // Pressure equation inner solve() loop max. iterations count
    std::size_t pressure_max_iter = 5;

    // Minimum residual for pressure equation
    double pressure_residual = 1e-7;
};


class PRIME : public IPressureLinked {
  public:
    PRIME(PRIMEParameters parameters);
    void step(std::span<eqn::Momentum*> momentum_predictors,
              field::Velocity& U,
              field::Velocity& mdot,
              field::Pressure& p) override;

  private:
    PRIMEParameters _params;
};

void solveMomentumExplicitly(std::span<eqn::Momentum*> momentum_predictors);

} // namespace prism::algo
