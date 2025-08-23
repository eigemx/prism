#pragma once

#include <cctype>
#include <span>

#include "algorithm.h"
#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"

namespace prism::algo {
struct SIMPLEParameters {
    // Under-relaxation factor for momentum predictor
    double momentum_urf = 0.7;

    // Momentum predictor inner solve() loop max. iterations count
    std::size_t momentum_max_iter = 3;

    // Minimum final solver residual for momentum predictor
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

class IncompressibleSIMPLE : public IPressureLinked {
  public:
    IncompressibleSIMPLE(SIMPLEParameters parameters);
    void step(std::span<eqn::Momentum*> momentum_predictors,
              SharedPtr<field::Velocity>& U,
              SharedPtr<field::Velocity>& mdot,
              SharedPtr<field::Pressure>& p) override;

  private:
    SIMPLEParameters _params;
};

void solveMomentumImplicitly(SIMPLEParameters params,
                             std::span<eqn::Momentum*> momentum_predictors);

void constrainPPrime(SharedPtr<field::Pressure>& pprime);

auto pressureEquationCoeffsTensor(std::span<eqn::Momentum*> momentum_predictors,
                                  const SharedPtr<field::Pressure>& p)
    -> SharedPtr<field::Tensor>;

auto solvePressureEquation(SIMPLEParameters params,
                           std::span<eqn::Momentum*> momentum_predictors,
                           SharedPtr<field::Velocity>& U,
                           SharedPtr<field::Velocity>& mdot,
                           SharedPtr<field::Pressure>& p)
    -> std::pair<SharedPtr<field::Pressure>, SharedPtr<field::Tensor>>;

void correctFields(SharedPtr<field::Velocity>& U,
                   SharedPtr<field::Velocity>& mdot,
                   SharedPtr<field::Pressure>& p,
                   const SharedPtr<field::Tensor>& D,
                   SharedPtr<field::Pressure>& pprime,
                   double pressure_urf);
} // namespace prism::algo
