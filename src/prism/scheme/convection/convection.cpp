#include "convection.h"
#include "prism/field/scalar.h"

namespace prism::scheme::convection {

IAppliedConvection::IAppliedConvection(const SharedPtr<field::IVector>& U,
                                       const SharedPtr<field::Scalar>& phi)
    : _U(U), IFullScheme(phi) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

void IAppliedConvection::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field()->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // since face is owned by `owner`, the flux will be in the intended way, as the face normal
    // points from cell center to face center, so no need to flup the sign of the flux.
    const f64 mdot_f = _U->fluxAtFace(face);

    auto [a_C, a_N, b] = interpolate(mdot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-mdot_f, neighbor, owner, face); // NOLINT

    this->insert(owner_id, owner_id, a_C);
    this->insert(owner_id, neighbor_id, a_N);

    this->insert(neighbor_id, neighbor_id, x_C);
    this->insert(neighbor_id, owner_id, x_N);

    this->rhs(owner_id) += b;
    this->rhs(neighbor_id) += s;
}

void IAppliedConvection::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::convection::IAppliedConvection",
                                           *this);
}

auto Upwind::interpolate(f64 m_dot,
                         const mesh::Cell& cell,     // NOLINT
                         const mesh::Cell& neighbor, // NOLINT
                         const mesh::Face& face)     // NOLINT
    -> detail::CoeffsTriplet {
    // when `cell` is the upwind cell
    f64 phi_tilde = phiTilde(field(), cell, neighbor);
    f64 phi_upwind = phiAtDummyUpwind(field(), cell, neighbor);
    const auto [l_plus, k_plus] = weightsNVF(phi_tilde);
    f64 a_C = std::max(m_dot, 0.0) * l_plus;
    f64 a_N = std::max(m_dot, 0.0) * k_plus;
    f64 b = -std::max(m_dot, 0.0) * (1 - l_plus - k_plus) * phi_upwind;

    // when `neighbor` is the upwind cell
    phi_tilde = phiTilde(field(), neighbor, cell);          // NOLINT
    phi_upwind = phiAtDummyUpwind(field(), neighbor, cell); // NOLINT
    const auto [l_minus, k_minus] = weightsNVF(phi_tilde);
    a_C += -std::max(-m_dot, 0.0) * k_minus;
    a_N += -std::max(-m_dot, 0.0) * l_minus;
    b += std::max(-m_dot, 0.0) * (1 - l_minus - k_minus) * phi_upwind;
    return {a_C, a_N, b};
}

auto Upwind::weightsNVF(f64 phi_tilde) const noexcept // NOLINT
    -> WeightsNVF {
    return {1.0, 0.0};
}

auto CentralDifference::weightsNVF(f64 phi_tilde) const noexcept // NOLINT
    -> Upwind::WeightsNVF {
    return {0.5, 0.5};
}

auto LinearUpwind::weightsNVF(f64 phi_tilde) const noexcept // NOLINT
    -> Upwind::WeightsNVF {
    return {1.5, 0.0};
}

auto QUICK::weightsNVF(f64 phi_tilde) const noexcept // NOLINT
    -> Upwind::WeightsNVF {
    return {0.75, 3.0 / 8.0};
}

auto FROMM::weightsNVF(f64 phi_tilde) const noexcept // NOLINT
    -> Upwind::WeightsNVF {
    return {1.0, 0.25};
}

auto MINMOD::weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < 0.5) {
        return {1.5, 0.0};
    }
    if (phi_tilde >= 0.5 && phi_tilde < 1.0) {
        return {0.5, 0.5};
    }
    return {1.0, 0};
}

auto MUSCL::weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < 0.25) {
        return {2, 0.0};
    }
    if (phi_tilde >= 0.25 && phi_tilde < 0.75) {
        return {1, 0.25};
    }
    if (phi_tilde >= 0.75 && phi_tilde < 1.0) {
        return {0, 1.0};
    }
    return {1.0, 0};
}

auto SMART::weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < (1 / 6)) {
        return {4, 0.0};
    }
    if (phi_tilde >= (1 / 6) && phi_tilde < (5 / 6)) {
        return {0.75, 3 / 8};
    }
    if (phi_tilde >= (5 / 6) && phi_tilde < 1.0) {
        return {0, 1.0};
    }
    return {1.0, 0};
}


} // namespace prism::scscheme::heme::convection
