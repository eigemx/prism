#pragma once

#include <cmath>

#include "boundary.h"
#include "convection_boundary.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/log.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/utilities.h"
#include "prism/operations/operations.h"
#include "prism/types.h"
#include "scheme.h"

namespace prism::scheme::convection {

// TODO: density should be a template parameter with a default value for field::UniformScalar
// TODO: we should have only one field::Vector allowed in IConvection schemes (rho * V) or just V

namespace detail {
// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double ownerCoeff {0.0};    // cell
    double neighborCoeff {0.0}; // neighbor
    double rhs {0.0};           // source
};
} // namespace detail

// Finite volume scheme for the discretization of the convection term
template <typename Field>
class IConvection : public IFullScheme<Field>,
                    public prism::boundary::BHManagersProvider<
                        boundary::ISchemeBoundaryHandler<IConvection<Field>>> {
  public:
    IConvection(field::Scalar rho, field::Velocity U, Field phi);

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

    auto inline field() -> Field override { return _phi; }
    auto inline U() -> const field::Velocity& { return _U; }
    auto inline rho() -> const field::Scalar& { return _rho; }


  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;

    void applyInterior(const mesh::Face& face) override;
    void applyBoundary();

    field::Scalar _rho;
    field::Velocity _U;
    Field _phi;
};

// Central difference scheme
template <typename F>
class CentralDifference : public IConvection<F> {
  public:
    CentralDifference(field::Scalar rho, field::Velocity U, F phi)
        : IConvection<F>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Upwind scheme
template <typename F>
class Upwind : public IConvection<F> {
  public:
    Upwind(field::Scalar rho, field::Velocity U, F phi) : IConvection<F>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Second order upwind scheme
template <typename F>
class SecondOrderUpwind : public IConvection<F> {
  public:
    SecondOrderUpwind(field::Scalar rho, field::Velocity U, F phi)
        : IConvection<F>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// QUICK scheme
template <typename F>
class QUICK : public IConvection<F> {
  public:
    QUICK(field::Scalar rho, field::Velocity U, F phi) : IConvection<F>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

template <typename Field>
IConvection<Field>::IConvection(field::Scalar rho, field::Velocity U, Field phi)
    : _rho(std::move(rho)),
      _U(std::move(U)),
      _phi(std::move(phi)),
      IFullScheme<Field>(phi.mesh().nCells()) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<scheme::boundary::Empty<Scheme>>();
    this->boundaryHandlersManager().template addHandler<scheme::boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<scheme::boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<scheme::boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<scheme::boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<scheme::boundary::NoSlip<Scheme>>();
}

template <typename Field>
void IConvection<Field>::apply() {
    applyBoundary();

    for (const auto& iface : _phi.mesh().interiorFaces()) {
        applyInterior(iface);
    }

    this->collect();
}

template <typename Field>
void IConvection<Field>::applyInterior(const mesh::Face& face) {
    const auto& mesh = _phi.mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    const Vector3d& S_f = mesh::outwardAreaVector(face, owner);
    const Vector3d U_f = _U.valueAtFace(face);
    const double rho_f = _rho.valueAtFace(face);
    const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

    auto [a_C, a_N, b] = interpolate(m_dot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-m_dot_f, neighbor, owner, face); // NOLINT

    this->insert(owner_id, owner_id, a_C);
    this->insert(owner_id, neighbor_id, a_N);

    this->insert(neighbor_id, neighbor_id, x_C);
    this->insert(neighbor_id, owner_id, x_N);

    this->rhs(owner_id) += b;
    this->rhs(neighbor_id) += s;
}

template <typename Field>
void IConvection<Field>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::convection::IConvection", *this);
}

template <typename F>
auto CentralDifference<F>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradScheme()->gradAtFace(face);
    const Vector3d d_Cf = face.center() - cell.center();
    const double a_C = std::max(m_dot, 0.0);
    double b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    const double a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

template <typename F>
auto Upwind<F>::interpolate(double m_dot,
                            const mesh::Cell& cell,     // NOLINT
                            const mesh::Cell& neighbor, // NOLINT
                            const mesh::Face& face)     // NOLINT
    -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const double a_C = std::max(m_dot, 0.0);
    // in case 'neighbor' is the upstream cell
    const double a_N = -std::max(-m_dot, 0.0);

    return {a_C, a_N, 0.0};
}

template <typename F>
auto SecondOrderUpwind<F>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradAtFace(face);
    const Vector3d cell_grad_phi = this->field().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->field().gradAtCell(neighbor);

    const Vector3d d_Cf = face.center() - cell.center();
    // auto correction = d_Cf.dot((2 * cell_grad_phi) - face_grad_phi);
    auto correction = cell_grad_phi.dot(d_Cf);

    const double a_C = std::max(m_dot, 0.0);
    const double b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    // correction = d_Nf.dot((2 * neighbor_grad_phi) - face_grad_phi);
    correction = neighbor_grad_phi.dot(d_Nf);

    const double a_N = -std::max(-m_dot, 0.0);
    const double b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}

template <typename F>
auto QUICK<F>::interpolate(double m_dot,
                           const mesh::Cell& cell,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->field().gradScheme()->gradAtFace(face);
    const Vector3d cell_grad_phi = this->field().gradScheme()->gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->field().gradScheme()->gradAtCell(neighbor);

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

namespace prism::scheme::boundary {
template <typename F>
void Fixed<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                              const mesh::BoundaryPatch& patch) {
    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const double phi_wall = patch.getScalarBoundaryCondition(phi.name());

        const Vector3d& S_f = face.area_vector();
        const Vector3d U_f = scheme.U().valueAtFace(face);

        // TODO: check if this is correct
        // TODO: warn user if mass flow rate is not entering the patch (negative flow rate)
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.

        // in case owner cell is an upstream cell
        // TODO: face value should be the same regardless of the upstream cell mass flow rate, so,
        // applying the right hand side should not depend on a std::max() call and should be
        // definite.
        const std::size_t cell_id = owner.id();
        scheme.insert(cell_id, cell_id, std::max(m_dot_f, 0.0));
        scheme.rhs(cell_id) += std::max(-m_dot_f * phi_wall, 0.0);
    }
}

template <typename F>
void NoSlip<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                               const mesh::BoundaryPatch& patch) {
    Fixed<convection::IConvection<F>> fixed;
    return fixed.apply(scheme, patch);
}

template <typename F>
void Outlet<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                               const mesh::BoundaryPatch& patch) {
    _n_reverse_flow_faces = 0;

    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const std::size_t cell_id = owner.id();

        // face area vector
        const Vector3d& S_f = face.area_vector();

        // use owner cell velocity as the velocity at the outlet face centroid
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

        if (m_dot_f < 0.0) {
            _n_reverse_flow_faces++;
        }
        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.

        // Because this is an outlet, owner `cell` is an upstream cell
        scheme.insert(cell_id, cell_id, m_dot_f);
    }

    if (_n_reverse_flow_faces > 0) {
        log::warn("Reverse flow detected in {} faces in outlet flow patch '{}'",
                  _n_reverse_flow_faces,
                  patch.name());
    }
}
} // namespace prism::scheme::boundary
