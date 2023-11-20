#include "relax.h"

#include <spdlog/spdlog.h>

#include "prism/constants.h"

namespace prism::solver {
void ImplicitUnderRelaxation::pre_relax(TransportEquation& eqn, double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        spdlog::debug(
            "ImplicitUnderRelaxation::pre_relax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& A = eqn.matrix();
    auto& b = eqn.rhs();
    const auto& phi = eqn.field().data();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= lambda;
    b += (1 - lambda) * A.diagonal().cwiseProduct(phi);
}

void ExplicitUnderRelaxation::post_relax(TransportEquation& eqn, double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        spdlog::debug(
            "ImplicitUnderRelaxation::pre_relax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& phi = eqn.field().data();
    const auto& phi_old = eqn.field_prev_iter().data();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda * (phi - phi_old));
}

} // namespace prism::solver