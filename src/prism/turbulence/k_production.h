#pragma once

#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::turbulence {

// Pk = tau_r : gradU
// where tau_r = mu_t * twoSymm(gradU) - (2/3) * rho * k * I
//
// For incompressible flow (div(U) = 0), the -(2/3)*rho*k*div(U) term vanishes,
// but is included for correctness and to keep the derivation explicit.
auto kProductionIncompressible(const SharedPtr<field::Scalar>& nu_t,
                               const SharedPtr<field::Scalar>& k,
                               const SharedPtr<field::Velocity>& U) -> SharedPtr<field::Scalar>;

// Compressible form adds the mu_t dilatation viscosity term to tau_r:
// tau_r = ... - (2/3) * mu_t * div(U) * I
// Pk   += -(2/3) * mu_t * div(U)^2
auto kProductionCompressible(const SharedPtr<field::Scalar>& mu_t,
                             const SharedPtr<field::Scalar>& rho,
                             const SharedPtr<field::Scalar>& k,
                             const SharedPtr<field::Velocity>& U) -> SharedPtr<field::Scalar>;

} // namespace prism::turbulence
