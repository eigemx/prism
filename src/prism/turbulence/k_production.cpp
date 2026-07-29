#include "k_production.h"

#include "prism/field/scalar.h"
#include "prism/operations/operations.h"

namespace prism::turbulence {

auto kProductionIncompressible(const SharedPtr<field::Scalar>& nu_t,
                               const SharedPtr<field::Scalar>& k,
                               const SharedPtr<field::Velocity>& U) -> SharedPtr<field::Scalar> {
    const auto& mesh = nu_t->mesh();
    const auto& nu_t_vals = nu_t->values();
    const auto& k_vals = k->values();
    VectorXd pk_vals(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        auto i = cell.id();
        Tensor3d gradU = U->gradAtCell(cell);
        Tensor3d tau_r = nu_t_vals[i] * ops::twoSymm(gradU);
        tau_r -= (2.0 / 3.0) * k_vals[i] * Matrix3d::Identity();
        pk_vals[i] = ops::doubleContraction(tau_r, gradU);
    }

    return makeShared<field::Scalar>("Pk", mesh, std::move(pk_vals));
}

auto kProductionCompressible(const SharedPtr<field::Scalar>& mu_t,
                             const SharedPtr<field::Scalar>& rho,
                             const SharedPtr<field::Scalar>& k,
                             const SharedPtr<field::Velocity>& U) -> SharedPtr<field::Scalar> {
    const auto& mesh = mu_t->mesh();
    const auto& mu_t_vals = mu_t->values();
    const auto& rho_vals = rho->values();
    const auto& k_vals = k->values();
    VectorXd pk_vals(mesh->cellCount());
    const auto& divU_vals = ops::div(*U);

    for (const auto& cell : mesh->cells()) {
        auto i = cell.id();
        Tensor3d gradU = U->gradAtCell(cell);
        Tensor3d I = Matrix3d::Identity();
        Tensor3d tau_r = mu_t_vals[i] * ops::twoSymm(gradU);
        tau_r -= (2.0 / 3.0) * rho_vals[i] * k_vals[i] * I;
        tau_r -= (2.0 / 3.0) * mu_t_vals[i] * divU_vals[i] * I;
        pk_vals[i] = ops::doubleContraction(tau_r, gradU);
    }

    return makeShared<field::Scalar>("Pk", mesh, std::move(pk_vals));
}

} // namespace prism::turbulence
