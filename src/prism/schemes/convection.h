#pragma once

#include <spdlog/spdlog.h>

#include <cmath>
#include <cstddef>

#include "fvscheme.h"
#include "prism/exceptions.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "prism/schemes/boundary.h"
#include "prism/types.h"


namespace prism::convection {
// forward declarations
template <typename GradScheme>
class IConvection;
} // namespace prism::convection

namespace prism::boundary {
template <typename G>
class Fixed<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename G>
class Empty<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme,
               const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "empty"; }
};

template <typename G>
class Symmetry<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme,
               const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename G>
class Outlet<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};
} // namespace prism::boundary

namespace prism::convection {

namespace detail {
// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {}; // cell
    double a_N {}; // neighbor
    double b {};   // source
};

auto inline face_mass_flow_rate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}
} // namespace detail

// Finite volume scheme for the discretization of the convection term
template <typename GradScheme = gradient::LeastSquares>
class IConvection : public FVScheme<field::Scalar> {
  public:
    IConvection(field::Scalar rho, field::Vector U, field::Scalar phi);
    void apply() override;
    auto inline field() -> std::optional<field::Scalar> override { return _phi; }

    using BCManager = boundary::BoundaryHandlersManager<IConvection<GradScheme>>;
    auto bc_manager() -> BCManager&;

  protected:
    auto inline grad_scheme() -> GradScheme& { return _gradient_scheme; }
    auto inline n_reverse_flow_faces() -> std::size_t& { return _n_reverse_flow_faces; }

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;

    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face);

    const field::Scalar _rho;
    const field::Vector _U;
    const field::Scalar _phi;
    GradScheme _gradient_scheme;

    std::size_t _n_reverse_flow_faces {0};
    BCManager _bc_manager;
};

// Central difference scheme
template <typename G = gradient::LeastSquares>
class CentralDifference : public IConvection<G> {
  public:
    CentralDifference(field::Scalar& rho, field::Vector& U, field::Scalar& phi)
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
    Upwind(field::Scalar rho, field::Vector U, field::Scalar phi) : IConvection<G>(rho, U, phi) {}

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
    SecondOrderUpwind(field::Scalar rho, field::Vector U, field::Scalar phi)
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
    QUICK(field::Scalar rho, field::Vector U, field::Scalar phi) : IConvection<G>(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

template <typename G>
IConvection<G>::IConvection(field::Scalar rho, field::Vector U, field::Scalar phi)
    : _rho(std::move(rho)),
      _U(std::move(U)),
      _phi(phi),
      _gradient_scheme(phi),
      FVScheme(phi.mesh().n_cells()) {}

template <typename G>
void IConvection<G>::apply() {
    _n_reverse_flow_faces = 0;

    for (const auto& bface : _phi.mesh().boundary_faces()) {
        apply_boundary(bface);
    }

    for (const auto& iface : _phi.mesh().interior_faces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    collect();

    if (_n_reverse_flow_faces > 0) {
        spdlog::warn(
            "convection::AbstractConvection::apply(): "
            "Reverse flow detected at {} faces in outlet boundary patch(es). "
            "This may cause the solution to diverge.",
            _n_reverse_flow_faces);
    }
}

template <typename G>
void IConvection<G>::apply_interior(const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const mesh::Cell& neighbor = _phi.mesh().cell(face.neighbor().value());

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    const Vector3d& S_f = mesh::outward_area_vector(face, owner);

    const Vector3d U_f = _U.value_at_face(face);
    const double rho_f = _rho.value_at_face(face);
    const double m_dot_f = detail::face_mass_flow_rate(rho_f, U_f, S_f);

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
void IConvection<G>::apply_boundary(const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const auto& boundary_patch = _phi.mesh().face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.kind()) {
        case mesh::BoundaryConditionKind::Empty:
        case mesh::BoundaryConditionKind::Symmetry: {
            return;
        }

        case mesh::BoundaryConditionKind::Fixed:
        case mesh::BoundaryConditionKind::VelocityInlet: {
            apply_boundary_fixed(owner, face);
            return;
        }

        case mesh::BoundaryConditionKind::Outlet: {
            apply_boundary_outlet(owner, face);
            return;
        }

            // TODO: Implement FixedGradient
        default:
            throw error::NonImplementedBoundaryCondition(
                "convection::AbstactConvection::apply_boundary()",
                boundary_patch.name(),
                boundary_condition.kind_string());
    }
}

template <typename G>
void IConvection<G>::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _phi.mesh().face_boundary_patch(face);
    const double phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    const Vector3d& S_f = face.area_vector();
    const Vector3d U_f = _U.value_at_face(face);

    // TODO: check if this is correct
    const double rho_f = _rho.value_at_cell(cell);
    const double m_dot_f = detail::face_mass_flow_rate(rho_f, U_f, S_f);

    // TODO: this assumes an upwind based scheme, this is wrong for central schemes
    // and should be generalized to work for all schemes.

    // in case owner cell is an upstream cell
    const std::size_t cell_id = cell.id();
    insert(cell_id, cell_id, std::max(m_dot_f, 0.0));
    rhs(cell_id) += std::max(-m_dot_f * phi_wall, 0.0);
}

template <typename G>
void IConvection<G>::apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face) {
    const std::size_t cell_id = cell.id();

    // face area vector
    const Vector3d& S_f = face.area_vector();

    // use owner cell velocity as the velocity at the outlet face centroid
    const Vector3d U_f = _U.value_at_face(face);
    const double rho_f = _rho.value_at_cell(cell);
    const double m_dot_f = detail::face_mass_flow_rate(rho_f, U_f, S_f);

    if (m_dot_f < 0.0) {
        n_reverse_flow_faces()++;
    }
    // TODO: this assumes an upwind based scheme, this is wrong for central schemes
    // and should be generalized to work for all schemes.

    // Because this is an outlet, owner `cell` is an upstream cell
    insert(cell_id, cell_id, m_dot_f);
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
    const Vector3d face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradient_at_cell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradient_at_cell(neighbor);

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
    const Vector3d face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradient_at_cell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradient_at_cell(neighbor);

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


} // namespace prism::convection
