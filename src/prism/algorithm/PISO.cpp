#include "PISO.h"

#include "PRIME.h"
#include "SIMPLE.h"

namespace prism::algo {

PISO::PISO(PISOControls controls) : _controls(controls) {}

auto PISO::step(std::span<eqn::Momentum*> momentum_predictors,
                SharedPtr<field::Velocity>& U,
                SharedPtr<field::Velocity>& mdot,
                SharedPtr<field::Pressure>& p) -> std::vector<report::Entry> {
    auto reports = std::vector<report::Entry>();

    SIMPLEControls simple_controls = {
        .momentum_urf = _controls.momentum_urf,
        .momentum_max_iter = _controls.momentum_max_iter,
        .momentum_residual = _controls.momentum_residual,
        .momentum_residual_stop = _controls.momentum_residual_stop,
        .non_ortho_correctors = _controls.non_ortho_correctors,
        .pressure_urf = _controls.pressure_urf,
        .pressure_max_iter = _controls.pressure_max_iter,
        .pressure_residual = _controls.pressure_residual,
        .p_ref_cell = _controls.p_ref_cell,
        .p_ref_value = _controls.p_ref_value,
    };

    PRIMEControls prime_controls = {
        .non_ortho_correctors = _controls.non_ortho_correctors,
        .pressure_urf = _controls.pressure_urf,
        .pressure_max_iter = _controls.pressure_max_iter,
        .pressure_residual = _controls.pressure_residual,
        .p_ref_cell = _controls.p_ref_cell,
        .p_ref_value = _controls.p_ref_value,
    };

    for (std::size_t i = 0; i < _controls.outer_iterations; ++i) {
        // Step 1: SIMPLE — implicit momentum predictor
        if (_controls.momentum_implicit_steps > 0) {
            for (std::size_t j = 0; j < _controls.momentum_implicit_steps; ++j) {
                auto simple_reports =
                    IncompressibleSIMPLE(simple_controls).step(momentum_predictors, U, mdot, p);
                reports.insert(reports.end(), simple_reports.begin(), simple_reports.end());
            }
        }

        // Steps 2-N: PRIME correctors (pressure_correction_steps includes the SIMPLE step)
        for (std::size_t j = 1; j < _controls.pressure_correction_steps; ++j) {
            auto prime_reports = PRIME(prime_controls).step(momentum_predictors, U, mdot, p);
            reports.insert(reports.end(), prime_reports.begin(), prime_reports.end());
        }
    }
    return reports;
}
} // namespace prism::algo
