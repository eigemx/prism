#include "relax.h"

#include "../print.h"

namespace prism::solver {
void ImplicitUnderRelaxation::pre_relax(Equation& eqn) const {
    auto& A = eqn.matrix();
    auto& b = eqn.rhs();

    const auto& phi = eqn.field().data();

    const auto lambda_ = lambda();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    // first apply the rhs update, before changing A.diagonal()
    b += ((1 - lambda_) / lambda_) * A.diagonal() * phi;
    // update diagonal coefficients
    A.diagonal() /= lambda_;
}

void ExplicitUnderRelaxation::post_relax(Equation& eqn) const {
    auto& phi = eqn.field().data();
    const auto& phi_old = eqn.field_prev_iter().data();
    const auto lambda_ = lambda();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda_ * (phi - phi_old));
}

} // namespace prism::solver