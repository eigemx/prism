#include "PRIME.h"

#include "SIMPLE.h"
#include "prism/equation/transport.h"

namespace prism::algo {

PRIME::PRIME(PRIMEParameters parameters) : _params(parameters) {}

void PRIME::step(std::span<eqn::Momentum*> momentum_predictors,
                      SharedPtr<field::Velocity>& U,
                      SharedPtr<field::Velocity>& mdot,
                      SharedPtr<field::Pressure>& p) {
    if (momentum_predictors.size() != 2 && momentum_predictors.size() != 3) {
        throw std::runtime_error(
            fmt::format("prism::algo::PRIME::step() expects 2 or 3 momentum predictors, not {}",
                        momentum_predictors.size()));
    }

    // solve momentum equations explicitly
    log::info("prism::algo::PRIME::step(): solving momentum equations explicitly");
    solveMomentumExplicitly(momentum_predictors);

    // solve pressure equation
    SIMPLEParameters SIMPLE_params {.pressure_urf = _params.pressure_urf,
                                    .pressure_max_iter = _params.pressure_max_iter,
                                    .pressure_residual = _params.pressure_residual};
    log::info("prism::algo::PRIME::step(): solving pressure equation");
    auto [pprime, D] = solvePressureEquation(SIMPLE_params, momentum_predictors, U, mdot, p);

    correctFields(U, mdot, p, D, pprime, _params.pressure_urf);
}

void solveMomentumExplicitly(std::span<eqn::Momentum*> momentum_predictors) {
    for (auto* eqn : momentum_predictors) {
        eqn->updateCoeffs();
        eqn->relax();

        auto U = eqn->field();
        const auto& A = eqn->matrix();
        const auto& b = eqn->rhs();
        const auto& ac = A.diagonal();

        const auto H = A * U->values() - ac.cwiseProduct(U->values());
        U->values() = (-H + b).cwiseQuotient(ac);

        eqn->zeroOutCoeffs();
    }
}
} // namespace prism::algo
