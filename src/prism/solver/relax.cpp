#include "relax.h"

#include "../print.h"

namespace prism::solver {
void ImplicitUnderRelaxation::pre_relax(Equation& eqn, double lambda) const {
    auto& A = eqn.matrix();
    auto& b = eqn.rhs();
    const auto& phi = eqn.field().data();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= lambda;
    b += (1 - lambda) * A.diagonal().cwiseProduct(phi);
}

void ExplicitUnderRelaxation::post_relax(Equation& eqn, double lambda) const {
    auto& phi = eqn.field().data();
    const auto& phi_old = eqn.field_prev_iter().data();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda * (phi - phi_old));
}

} // namespace prism::solver