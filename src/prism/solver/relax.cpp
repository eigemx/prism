#include "relax.h"

namespace prism::solver {
void ImplicitUnderRelaxation::relax(Equation& eqn) const {
    auto& A = eqn.coeff_matrix();
    auto& b = eqn.rhs_vector();
    const auto& phi_old = eqn.scalar_field_old().data();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= _lambda;
    b += ((1 - _lambda) / _lambda) * A.diagonal() * phi_old;
}
} // namespace prism::solver