#pragma once

#include <cctype>
#include <span>

#include "algorithm.h"
#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"

namespace prism::algo {
struct SIMPLEParameters {
    // Number of non-orthogonal correctors for pressure equation
    std::size_t non_ortho_correctors = 2;

    // Under-relaxation factor for momentum predictor
    double momentum_urf = 0.7;

    // Momentum predictor inner solve() loop max. iterations count
    std::size_t momentum_max_iter = 3;

    // Minimum residual for momentum predictor
    double momentum_residual = 1e-7;

    // Under-relaxation factor for pressure field
    double pressure_urf = 0.3;

    // Pressure equation inner solve() loop max. iterations count
    std::size_t pressure_max_iter = 5;

    // Minimum residual for pressure equation
    double pressure_residual = 1e-7;
};

class IncompressibleSIMPLE : public IPressureLinked {
  public:
    IncompressibleSIMPLE(SIMPLEParameters parameters);
    void step(std::span<eqn::Momentum*> momentum_predictors,
              field::Velocity& U,
              field::Velocity& mdot,
              field::Pressure& p) override;

  private:
    SIMPLEParameters _params;
};
void constrainPPrime(field::Pressure& pprime);

auto pressureEquationCoeffsTensor(std::span<eqn::Momentum*> momentum_predictors,
                                  const field::Pressure& p) -> field::Tensor;

auto solvePressureEquation(SIMPLEParameters params,
                           std::span<eqn::Momentum*> momentum_predictors,
                           field::Velocity& U,
                           field::Velocity& mdot,
                           const field::Pressure& p) -> std::pair<field::Pressure, field::Tensor>;

void correctFields(field::Velocity& U,
                   field::Velocity& mdot,
                   field::Pressure& p,
                   const field::Tensor& D,
                   const field::Pressure& pprime,
                   double pressure_urf);
} // namespace prism::algo
