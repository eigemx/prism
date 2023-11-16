#pragma once

#include <cmath>

#include "fvscheme.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"


namespace prism::convection {
// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {}; // cell
    double a_N {}; // neighbor
    double b {};   // source
};

auto inline face_mass_flow_rate(double rho, const prism::Vector3d& U, const prism::Vector3d& S)
    -> double {
    return rho * U.dot(S);
}

// Finite volume scheme for the discretization of the convection term
template <typename GradScheme>
class AbstractConvection : public FVScheme {
  public:
    AbstractConvection(ScalarField& rho, VectorField& U, ScalarField& phi);

    void apply() override;

    auto inline field() -> std::optional<ScalarField> override { return _phi; }
    auto inline grad_scheme() -> GradScheme& { return _gradient_scheme; }

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> CoeffsTriplet = 0;

    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face);
    auto boundary_face_velocity(const mesh::Face& face) const -> Vector3d;

    ScalarField _rho;
    VectorField& _U;
    ScalarField _phi;
    const mesh::PMesh& _mesh;
    GradScheme _gradient_scheme;
};

// Central difference scheme
template <typename G = gradient::LeastSquares>
class CentralDifference : public AbstractConvection<G> {
  public:
    CentralDifference(ScalarField& rho, VectorField& U, ScalarField& phi)
        : AbstractConvection<G>(rho, U, phi) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> CoeffsTriplet override;
};


// Upwind scheme
template <typename G = gradient::LeastSquares>
class Upwind : public AbstractConvection<G> {
  public:
    Upwind(ScalarField& rho, VectorField& U, ScalarField& phi)
        : AbstractConvection<G>(rho, U, phi) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> CoeffsTriplet override;
};


// Second order upwind scheme
template <typename G = gradient::LeastSquares>
class SecondOrderUpwind : public AbstractConvection<G> {
  public:
    SecondOrderUpwind(ScalarField& rho, VectorField& U, ScalarField& phi)
        : AbstractConvection<G>(rho, U, phi) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> CoeffsTriplet override;
};


// QUICK scheme
template <typename G = gradient::LeastSquares>
class QUICK : public AbstractConvection<G> {
  public:
    QUICK(ScalarField& rho, VectorField& U, ScalarField& phi)
        : AbstractConvection<G>(rho, U, phi) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> CoeffsTriplet override;
};

template <typename G>
AbstractConvection<G>::AbstractConvection(ScalarField& rho, VectorField& U, ScalarField& phi)
    : _rho(rho),
      _U(U),
      _phi(phi),
      _mesh(phi.mesh()),
      _gradient_scheme(phi),
      FVScheme(phi.mesh().n_cells()) {}

template <typename G>
void AbstractConvection<G>::apply() {
    // TODO: this is repeated in all FVSchemes, we should move it to the base class
    for (const auto& bface : _mesh.boundary_faces()) {
        apply_boundary(_mesh.cell(bface.owner()), bface);
    }

    for (const auto& iface : _mesh.interior_faces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    collect();
}

template <typename G>
void AbstractConvection<G>::apply_interior(const mesh::Face& face) {
    const auto& owner = _mesh.cell(face.owner());
    const auto& neighbor = _mesh.cell(face.neighbor().value());

    auto owner_id = owner.id();
    auto neighbor_id = neighbor.id();

    auto S_f = mesh::outward_area_vector(face, owner);

    // interpolated velocity vector and density at face centroid
    auto g_c = mesh::geo_weight(owner, neighbor, face);

    // TODO: convert this to a utility function
    auto U_f = (g_c * _U[owner_id]) + ((1 - g_c) * _U[neighbor_id]);
    auto rho_f = (g_c * _rho[owner_id]) + ((1 - g_c) * _rho[neighbor_id]);

    auto m_dot_f = face_mass_flow_rate(rho_f, U_f, S_f);

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
void AbstractConvection<G>::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.bc_type()) {
        case mesh::BoundaryConditionType::Empty:
        case mesh::BoundaryConditionType::Symmetry: {
            return;
        }

        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            apply_boundary_fixed(cell, face);
            return;
        }

        case mesh::BoundaryConditionType::Outlet: {
            apply_boundary_outlet(cell, face);
            return;
        }

        default:
            throw std::runtime_error(fmt::format(
                "convection::ConvectionBase::apply_boundary(): "
                "Non-implemented boundary type for boundary patch: '{}' for field '{}'",
                boundary_patch.name(),
                _U.name()));
    }
}

template <typename G>
void AbstractConvection<G>::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    const auto& S_f = face.area_vector();
    const auto& U_f = boundary_face_velocity(face);

    // TODO: check if this is correct
    auto rho_f = _rho[face.owner()];
    auto m_dot_f = face_mass_flow_rate(rho_f, U_f, S_f);

    // TODO: this assumes an upwind based scheme, this is wrong for central schemes
    // and should be generalized to work for all schemes.

    // in case owner cell is an upstream cell
    insert(cell.id(), cell.id(), std::max(m_dot_f, 0.0));
    rhs(cell.id()) += std::max(-m_dot_f * phi_wall, 0.0);
}

template <typename G>
void AbstractConvection<G>::apply_boundary_outlet(const mesh::Cell& cell,
                                                  const mesh::Face& face) {
    auto cell_id = cell.id();

    // face area vector
    const auto& S_f = face.area_vector();

    // use owner cell velocity as the velocity at the outlet face centroid
    const auto& U_f = boundary_face_velocity(face);

    auto rho_f = _rho[face.owner()];

    auto m_dot_f = face_mass_flow_rate(rho_f, U_f, S_f);

    if (m_dot_f <= 0.0) {
        warn(
            fmt::format("convection::ConvectionBase::apply_boundary_outlet(): "
                        "Reverse flow detected at outlet boundary patch '{}'. "
                        "This may cause the solution to diverge.",
                        _mesh.face_boundary_patch(face).name()));
    }
    // TODO: this assumes an upwind based scheme, this is wrong for central schemes
    // and should be generalized to work for all schemes.

    // Because this is an outlet, owner `cell` is an upstream cell
    insert(cell_id, cell_id, m_dot_f);
}

template <typename G>
auto AbstractConvection<G>::boundary_face_velocity(const mesh::Face& face) const -> Vector3d {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_U.name());

    switch (boundary_condition.bc_type()) {
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            return boundary_patch.get_vector_bc(_U.name());
        }

        case mesh::BoundaryConditionType::Outlet:
        case mesh::BoundaryConditionType::Symmetry: {
            return _U[face.owner()];
        }

        default:
            throw std::runtime_error(fmt::format(
                "convection::IConvection::boundary_face_velocity(): "
                "Non-implemented boundary type for boundary patch: '{}' for field '{}'",
                boundary_patch.name(),
                _U.name()));
    }
}

template <typename G>
auto CentralDifference<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    auto face_grad_phi = this->grad_scheme().gradient_at_face(face);
    auto d_Cf = face.center() - cell.center();
    auto a_C = std::max(m_dot, 0.0);
    auto b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    auto a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

template <typename G>
auto Upwind<G>::interpolate(double m_dot,
                            const mesh::Cell& cell,     // NOLINT
                            const mesh::Cell& neighbor, // NOLINT
                            const mesh::Face& face)     // NOLINT
    -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    auto a_C = std::max(m_dot, 0.0);
    // in case 'neighbor' is the upstream cell
    auto a_N = -std::max(-m_dot, 0.0);

    auto b = 0.0;

    return {a_C, a_N, b};
}

template <typename G>
auto SecondOrderUpwind<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    const auto& face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const auto& cell_grad_phi = this->grad_scheme().gradient_at_cell(cell);
    const auto& neighbor_grad_phi = this->grad_scheme().gradient_at_cell(neighbor);

    auto d_Cf = face.center() - cell.center();
    //auto correction = d_Cf.dot((2 * cell_grad_phi) - face_grad_phi);
    auto correction = cell_grad_phi.dot(d_Cf);

    auto a_C = std::max(m_dot, 0.0);
    auto b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    //correction = d_Nf.dot((2 * neighbor_grad_phi) - face_grad_phi);
    correction = neighbor_grad_phi.dot(d_Nf);

    auto a_N = -std::max(-m_dot, 0.0);
    auto b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}

template <typename G>
auto QUICK<G>::interpolate(double m_dot,
                           const mesh::Cell& cell,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    const auto& face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const auto& cell_grad_phi = this->grad_scheme().gradient_at_cell(cell);
    const auto& neighbor_grad_phi = this->grad_scheme().gradient_at_cell(neighbor);

    auto d_Cf = face.center() - cell.center();
    auto correction = 0.5 * d_Cf.dot(cell_grad_phi + face_grad_phi);

    auto a_C = std::max(m_dot, 0.0);
    auto b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    correction = 0.5 * d_Nf.dot(neighbor_grad_phi + face_grad_phi);

    auto a_N = -std::max(-m_dot, 0.0);
    auto b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}


} // namespace prism::convection
