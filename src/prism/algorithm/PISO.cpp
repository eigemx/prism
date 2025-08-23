#include "PISO.h"

#include "PRIME.h"
#include "SIMPLE.h"

namespace prism::algo {

PISO::PISO(PISOParameters parameters) : _params(parameters) {}

void PISO::step(std::span<eqn::Momentum*> momentum_predictors,
                SharedPtr<field::Velocity>& U,
                SharedPtr<field::Velocity>& mdot,
                SharedPtr<field::Pressure>& p) {
    for (std::size_t i = 0; i < _params.pressure_correction_steps; ++i) {
        if (_params.momentum_implicit_steps > 0) {
            SIMPLEParameters simple_params = {.momentum_urf = _params.momentum_urf,
                                              .momentum_max_iter = _params.momentum_max_iter,
                                              .momentum_residual = _params.momentum_residual};

            for (std::size_t j = 0; j < _params.momentum_implicit_steps; ++j) {
                IncompressibleSIMPLE(simple_params).step(momentum_predictors, U, mdot, p);
            }
        }
        for (std::size_t j = 0; j < _params.pressure_correction_steps; ++j) {
            PRIMEParameters prime_params = {
                .non_ortho_correctors = _params.non_ortho_correctors,
                .pressure_urf = _params.pressure_urf,
                .pressure_max_iter = _params.pressure_max_iter,
                .pressure_residual = _params.pressure_residual,
            };
            PRIME(prime_params).step(momentum_predictors, U, mdot, p);
        }
    }
}
} // namespace prism::algo
