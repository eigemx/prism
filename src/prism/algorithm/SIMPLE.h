#pragma once

#include <span>

#include "algorithm.h"
#include "prism/equation/boundary.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/field/velocity.h"
#include "prism/report.h"
#include "prism/types.h"

namespace prism::algo {
struct SIMPLEParameters {
    // Under-relaxation factor for momentum predictor
    f64 momentum_urf = 0.7;

    // Momentum predictor inner solve() loop max. iterations count
    size_t momentum_max_iter = 3;

    // Minimum final solver residual for momentum predictor
    f64 momentum_residual = 1e-7;

    // Initial residual stopping criterion for momentum predictor
    f64 momentum_residual_stop = 1e-5;

    // Number of non-orthogonal correctors for pressure equation
    size_t non_ortho_correctors = 2;

    // Under-relaxation factor for pressure field
    f64 pressure_urf = 0.3;

    // Pressure equation inner solve() loop max. iterations count
    size_t pressure_max_iter = 5;

    // Minimum residual for pressure equation
    f64 pressure_residual = 1e-7;

    // Reference cell for singular (all-Neumann) pressure systems.
    // When set, the diagonal at this cell is doubled to break the null space.
    Optional<size_t> p_ref_cell = NullOption;
    f64 p_ref_value = 0.0;
};

class IncompressibleSIMPLE : public IPressureLinked {
  public:
    IncompressibleSIMPLE(SIMPLEParameters parameters);
    auto step(std::span<eqn::Momentum*> momentum_predictors,
              SharedPtr<field::Velocity>& U,
              SharedPtr<field::Velocity>& mdot,
              SharedPtr<field::Pressure>& p) -> std::vector<report::Entry> override;

  private:
    SIMPLEParameters _params;
};

auto solveImplicitMomentum(SIMPLEParameters params, std::span<eqn::Momentum*> momentum_predictors)
    -> std::vector<report::Entry>;

// Returns true if all pressure BCs are Neumann (closed domain requiring a reference cell).
auto needsPressureReference(const SharedPtr<field::Pressure>& pprime) -> bool;

void constrainPPrime(SharedPtr<field::Pressure>& pprime);

auto pressureEquationCoeffsTensor(std::span<eqn::Momentum*> momentum_predictors,
                                  const SharedPtr<field::Pressure>& p) -> SharedPtr<field::Tensor>;

struct PressureEquationResult {
    SharedPtr<field::Pressure> pprime;
    SharedPtr<field::Tensor> D;
    std::vector<report::Entry> reports;
};

auto solvePressureEquation(SIMPLEParameters params,
                           std::span<eqn::Momentum*> momentum_predictors,
                           SharedPtr<field::Velocity>& U,
                           SharedPtr<field::Velocity>& mdot,
                           SharedPtr<field::Pressure>& p) -> PressureEquationResult;

void correctFields(SharedPtr<field::Velocity>& U,
                   SharedPtr<field::Velocity>& mdot,
                   SharedPtr<field::Pressure>& p,
                   const SharedPtr<field::Tensor>& D,
                   SharedPtr<field::Pressure>& pprime,
                   double pressure_urf);
} // namespace prism::algo
