#include "PRIME.h"

#include "SIMPLE.h"
#include "prism/equation/transport.h"
#include "prism/log.h"

namespace prism::algo {

PRIME::PRIME(PRIMEControls controls) : _controls(controls) {}

auto PRIME::step(std::span<eqn::Momentum*> momentum_predictors,
                 SharedPtr<field::Velocity>& U,
                 SharedPtr<field::Velocity>& mdot,
                 SharedPtr<field::Pressure>& p) -> std::vector<report::Entry> {
    if (momentum_predictors.size() != 2 && momentum_predictors.size() != 3) {
        throw std::runtime_error(
            fmt::format("prism::algo::PRIME::step() expects 2 or 3 momentum predictors, not {}",
                        momentum_predictors.size()));
    }

    // solve momentum equations explicitly
    log::debug("prism::algo::PRIME::step(): solving momentum equations explicitly");
    solveExplicitMomentum(momentum_predictors);

    // solve pressure equation
    SIMPLEControls SIMPLE_controls {.pressure_urf = _controls.pressure_urf,
                                    .pressure_max_iter = _controls.pressure_max_iter,
                                    .pressure_residual = _controls.pressure_residual,
                                    .p_ref_cell = _controls.p_ref_cell,
                                    .p_ref_value = _controls.p_ref_value};
    log::debug("prism::algo::PRIME::step(): solving pressure equation");
    auto result = solvePressureEquation(SIMPLE_controls, momentum_predictors, U, mdot, p);

    correctFields(U, mdot, p, result.D, result.pprime, _controls.pressure_urf);
    return result.reports;
}

void solveExplicitMomentum(std::span<eqn::Momentum*> momentum_predictors) {
    for (auto* eqn : momentum_predictors) {
        eqn->updateCoeffs();
        eqn->relax();

        auto U = eqn->field();
        const auto& A = eqn->matrix();
        const auto& b = eqn->rhs();
        const auto& ac = A.diagonal();

        const auto H = A * U->values() - ac.cwiseProduct(U->values());
        U->values() = (-H + b).cwiseQuotient(ac);

        // Keep the relaxed matrix: pressureEquationCoeffsTensor() reuses its diagonal to build D.
    }
}
} // namespace prism::algo
