#pragma once

#include <spdlog/spdlog.h>

#include <cmath>

#include "convection_boundary.h"
#include "fvscheme.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "prism/schemes/boundary.h"
#include "prism/types.h"

namespace prism::scheme::convection {

namespace detail {
// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {0.0}; // cell
    double a_N {0.0}; // neighbor
    double b {0.0};   // source
};
} // namespace detail

// Finite volume scheme for the discretization of the convection term
template <typename GradScheme = gradient::LeastSquares>
class IConvection : public FVScheme<field::Scalar> {
  public:
    IConvection(field::Scalar rho, field::Velocity U, field::Scalar phi);

    void apply() override;

    auto inline field() -> std::optional<field::Scalar> override { return _phi; }
    auto inline U() -> const field::Velocity& { return _U; }
    auto inline rho() -> const field::Scalar& { return _rho; }

    using GradSchemeType = GradScheme;

    using BoundaryHandlersManager = prism::boundary::BoundaryHandlersManager<
        IConvection<GradScheme>,
        boundary::FVSchemeBoundaryHandler<IConvection<GradScheme>>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }

  protected:
    auto inline grad_scheme() -> GradScheme& { return _gradient_scheme; }

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;

    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    field::Scalar _rho;
    field::Velocity _U;
    field::Scalar _phi;
    GradScheme _gradient_scheme;
    BoundaryHandlersManager _bc_manager;
};

// Central difference scheme
template <typename G = gradient::LeastSquares>
class CentralDifference : public IConvection<G> {
  public:
    CentralDifference(field::Scalar& rho, field::Velocity& U, field::Scalar& phi)
        : IConvection<G>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Upwind scheme
template <typename G = gradient::LeastSquares>
class Upwind : public IConvection<G> {
  public:
    Upwind(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection<G>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Second order upwind scheme
template <typename G = gradient::LeastSquares>
class SecondOrderUpwind : public IConvection<G> {
  public:
    SecondOrderUpwind(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection<G>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// QUICK scheme
template <typename G = gradient::LeastSquares>
class QUICK : public IConvection<G> {
  public:
    QUICK(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection<G>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

template <typename G>
IConvection<G>::IConvection(field::Scalar rho, field::Velocity U, field::Scalar phi)
    : _rho(std::move(rho)),
      _U(std::move(U)),
      _phi(phi),
      _gradient_scheme(phi),
      FVScheme(phi.mesh().nCells()) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Symmetry<Scheme>>();
}

template <typename G>
void IConvection<G>::apply() {
    apply_boundary();

    for (const auto& iface : _phi.mesh().interiorFaces()) {
        apply_interior(iface);
    }

    collect();
}

template <typename G>
void IConvection<G>::apply_interior(const mesh::Face& face) {
    const auto& mesh = _phi.mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    const Vector3d& S_f = mesh::outward_area_vector(face, owner);
    const Vector3d U_f = _U.valueAtFace(face);
    const double rho_f = _rho.valueAtFace(face);
    const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

    auto [a_C, a_N, b] = interpolate(m_dot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-m_dot_f, neighbor, owner, face); // NOLINT

    insert(owner_id, owner_id, a_C);
    insert(owner_id, neighbor_id, a_N);

    insert(neighbor_id, neighbor_id, x_C);
    insert(neighbor_id, owner_id, x_N);

    rhs(owner_id) += b;
    rhs(neighbor_id) += s;
}

template <typename G>
void IConvection<G>::apply_boundary() {
    boundary::detail::apply_boundary("IConvection", *this);
}

template <typename G>
auto CentralDifference<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const Vector3d d_Cf = face.center() - cell.center();
    const double a_C = std::max(m_dot, 0.0);
    double b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    const double a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

template <typename G>
auto Upwind<G>::interpolate(double m_dot,
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

template <typename G>
auto SecondOrderUpwind<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradAtFace(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradAtCell(neighbor);

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

template <typename G>
auto QUICK<G>::interpolate(double m_dot,
                           const mesh::Cell& cell,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradAtFace(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradAtCell(neighbor);

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
