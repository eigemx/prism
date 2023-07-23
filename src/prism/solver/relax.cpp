#include "relax.h"

namespace prism::solver {
void ImplicitUnderRelaxation::pre_relax(Equation& eqn) const {
    auto& A = eqn.coeff_matrix();
    auto& b = eqn.rhs_vector();
    const auto& phi_old = eqn.scalar_field_old().data();
    const auto lambda_ = lambda();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= lambda_;
    b += ((1 - lambda_) / lambda_) * A.diagonal() * phi_old;
}

void ExplicitUnderRelaxation::post_relax(Equation& eqn) const {
    auto& phi = eqn.scalar_field().data();
    const auto& phi_old = eqn.scalar_field_old().data();
    const auto lambda_ = lambda();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda_ * (phi - phi_old));
}

} // namespace prism::solver