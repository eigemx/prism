#pragma once

#include <cctype>
#include <span>

#include "algorithm.h"

namespace prism::algo {

struct PISOParameters {
    // PISO outer loop iterations count
    std::size_t outer_iterations = 1;

    // Number of SIMPLE steps
    std::size_t momentum_implicit_steps = 1;

    // Number of pressure correction steps (PRIME steps).
    std::size_t pressure_correction_steps = 2;

    // Under-relaxation factor for momentum predictor
    double momentum_urf = 0.7;

    // Momentum predictor inner solve() loop max. iterations count
    std::size_t momentum_max_iter = 3;

    // Minimum residual for momentum predictor
    double momentum_residual = 1e-7;

    // Initial residual stopping criterion for momentum predictor
    double momentum_residual_stop = 1e-5;

    // Number of non-orthogonal correctors for pressure equation
    std::size_t non_ortho_correctors = 2;

    // Under-relaxation factor for pressure field
    double pressure_urf = 0.3;

    // Pressure equation inner solve() loop max. iterations count
    std::size_t pressure_max_iter = 5;

    // Minimum residual for pressure equation
    double pressure_residual = 1e-7;
};

class PISO : public IPressureLinked {
  public:
    PISO(PISOParameters parameters);
    void step(std::span<eqn::Momentum*> momentum_predictors,
              field::Velocity& U,
              field::Velocity& mdot,
              field::Pressure& p) override;

  private:
    PISOParameters _params;
};

} // namespace prism::algo
