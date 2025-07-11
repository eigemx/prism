#pragma once

#include <cmath>

#include "boundary.h"
#include "convection_boundary.h"
#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"
#include "scheme.h"

namespace prism::scheme::convection {

namespace detail {
// coefficients for the discretized convection equation
struct CoeffsTriplet {
    double ownerCoeff {0.0};    // cell
    double neighborCoeff {0.0}; // neighbor
    double rhs {0.0};           // source
};
} // namespace detail

// Basic base class for all convection schemes, without templating clutter.
class IConvection {};

// Finite volume scheme for the discretization of the convection term
template <field::IVectorBased ConvectiveField, typename Field>
class IAppliedConvection
    : public IConvection,
      public IFullScheme<Field>,
      public prism::boundary::BHManagerProvider<
          boundary::ISchemeBoundaryHandler<IAppliedConvection<ConvectiveField, Field>>> {
  public:
    IAppliedConvection(ConvectiveField U, Field phi);

    auto needsCorrection() const noexcept -> bool override { return true; }
    auto inline field() -> Field override { return _phi; }
    auto inline U() -> const ConvectiveField& { return _U; }

    using ConvectiveFieldType = ConvectiveField;
    using FieldType = Field;

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;
    ConvectiveField _U;
    Field _phi;

    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

// Concept for diffusion schemes that are based on IAppliedConvection.
template <typename T>
concept IAppliedConvectionBased =
    std::derived_from<T,
                      IAppliedConvection<typename T::ConvectiveFieldType, typename T::FieldType>>;

// Central difference scheme
template <field::IVectorBased ConvectiveField, typename Field>
class CentralDifference : public IAppliedConvection<ConvectiveField, Field> {
  public:
    CentralDifference(ConvectiveField U, Field phi)
        : IAppliedConvection<ConvectiveField, Field>(U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

// Upwind scheme
template <field::IVectorBased ConvectiveField, typename Field>
class Upwind : public IAppliedConvection<ConvectiveField, Field> {
  public:
    Upwind(ConvectiveField U, Field phi) : IAppliedConvection<ConvectiveField, Field>(U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

// Second order upwind (linear upwind) scheme
template <field::IVectorBased ConvectiveField, typename Field>
class LinearUpwind : public IAppliedConvection<ConvectiveField, Field> {
  public:
    LinearUpwind(ConvectiveField U, Field phi)
        : IAppliedConvection<ConvectiveField, Field>(U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

// QUICK scheme
template <field::IVectorBased ConvectiveField, typename Field>
class QUICK : public IAppliedConvection<ConvectiveField, Field> {
  public:
    QUICK(ConvectiveField U, Field phi) : IAppliedConvection<ConvectiveField, Field>(U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

template <field::IVectorBased ConvectiveField, typename Field>
IAppliedConvection<ConvectiveField, Field>::IAppliedConvection(ConvectiveField U, Field phi)
    /// TODO: check why _phi.mesh()->cellCount() is not working, as phi should be in a moved
    /// state. Also, we need to avoid std::move and just make the constructor take a const
    /// reference.
    : _U(std::move(U)), _phi(std::move(phi)), IFullScheme<Field>(phi.mesh()->cellCount()) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <field::IVectorBased ConvectiveField, typename Field>
void IAppliedConvection<ConvectiveField, Field>::applyInterior(const mesh::Face& face) {
    const auto& mesh = _phi.mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // since face is owned by `owner`, the flux will be in the intended way, as the face normal
    // points from cell center to face center, so no need to flup the sign of the flux.
    const double mdot_f = _U.fluxAtFace(face);

    auto [a_C, a_N, b] = interpolate(mdot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-mdot_f, neighbor, owner, face); // NOLINT

    this->insert(owner_id, owner_id, a_C);
    this->insert(owner_id, neighbor_id, a_N);

    this->insert(neighbor_id, neighbor_id, x_C);
    this->insert(neighbor_id, owner_id, x_N);

    this->rhs(owner_id) += b;
    this->rhs(neighbor_id) += s;
}

template <field::IVectorBased ConvectiveField, typename Field>
void IAppliedConvection<ConvectiveField, Field>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::convection::IAppliedConvection",
                                           *this);
}

template <field::IVectorBased ConvectiveField, typename Field>
auto CentralDifference<ConvectiveField, Field>::interpolate(double m_dot,
                                                            const mesh::Cell& cell,
                                                            const mesh::Cell& neighbor,
                                                            const mesh::Face& face)
    -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradAtFace(face);
    const Vector3d d_Cf = face.center() - cell.center();
    const double a_C = std::max(m_dot, 0.0);
    double b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    const double a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto Upwind<ConvectiveField, Field>::interpolate(double m_dot,
                                                 const mesh::Cell& cell,     // NOLINT
                                                 const mesh::Cell& neighbor, // NOLINT
                                                 const mesh::Face& face)     // NOLINT
    -> detail::CoeffsTriplet {
    const double a_C = std::max(m_dot, 0.0);
    const double a_N = -std::max(-m_dot, 0.0);
    return {a_C, a_N, 0.0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto LinearUpwind<ConvectiveField, Field>::interpolate(double m_dot,
                                                       const mesh::Cell& cell,
                                                       const mesh::Cell& neighbor,
                                                       const mesh::Face& face)
    -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradAtFace(face);
    const Vector3d cell_grad_phi = this->field().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->field().gradAtCell(neighbor);

    const Vector3d d_Cf = face.center() - cell.center();
    auto correction = d_Cf.dot((2 * cell_grad_phi) - face_grad_phi);
    // auto correction = cell_grad_phi.dot(d_Cf);

    const double a_C = std::max(m_dot, 0.0);
    const double b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    correction = d_Nf.dot((2 * neighbor_grad_phi) - face_grad_phi);
    // correction = neighbor_grad_phi.dot(d_Nf);

    const double a_N = -std::max(-m_dot, 0.0);
    const double b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto QUICK<ConvectiveField, Field>::interpolate(double m_dot,
                                                const mesh::Cell& cell,
                                                const mesh::Cell& neighbor,
                                                const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradAtFace(face);
    const Vector3d cell_grad_phi = this->field().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->field().gradAtCell(neighbor);

    const Vector3d d_Cf = face.center() - cell.center();
    auto correction = 0.5 * d_Cf.dot(cell_grad_phi + face_grad_phi);

    const double a_C = std::max(m_dot, 0.0);
    const double b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    correction = 0.5 * d_Nf.dot(neighbor_grad_phi + face_grad_phi);

    const double a_N = -std::max(-m_dot, 0.0);
    const double b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}
} // namespace prism::scheme::convection
