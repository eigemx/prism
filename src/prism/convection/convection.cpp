#include "convection.h"

#include "../mesh/utilities.h"

namespace prism::convection {
void ConvectionBase::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();

    // get adjacent cell sharing `face` with `cell`
    auto is_owned = face.is_owned_by(cell_id);
    auto adjacent_cell_id = is_owned ? face.neighbor().value() : face.owner();

    const auto& adj_cell = _mesh.cell(adjacent_cell_id);

    // density at cell centroid
    auto rho_c = _rho[cell_id];

    auto S_f = mesh::outward_area_vector(face, cell);

    // interpolated velocity vector at face centroid
    auto g_c = mesh::geo_weight(cell, adj_cell, face);
    const auto& U_f = (g_c * _U[cell_id]) + ((1 - g_c) * _U[adjacent_cell_id]);

    auto m_dot_f = face_mass_flow_rate(rho_c, U_f, S_f);

    auto [a_C, a_N, b] = interpolate(m_dot_f, cell, adj_cell, face, _gradient_scheme);


    // TODO: This assumes that velocity field is constant, this wrong because
    // the velocity field when solving the momentum equation will be different in
    // each iteration. This should be generalized to work for all schemes.
    // a possible workaround is to use zero out the coefficients matrix every
    // iteration in finalaize(), same goes for below member functions.
    matrix(cell_id, cell_id) += a_C;
    matrix(cell_id, adjacent_cell_id) += a_N;
    rhs(cell_id) += b;
}

void ConvectionBase::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
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

void ConvectionBase::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
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
    matrix(cell.id(), cell.id()) += std::max(m_dot_f, 0.0);

    rhs(cell.id()) += std::max(-m_dot_f * phi_wall, 0.0);
}

void ConvectionBase::apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face) {
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
    matrix(cell_id, cell_id) += m_dot_f;
}

auto ConvectionBase::boundary_face_velocity(const mesh::Face& face) const -> Vector3d {
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
                "convection::ConvectionBase::boundary_face_velocity(): "
                "Non-implemented boundary type for boundary patch: '{}' for field '{}'",
                boundary_patch.name(),
                _U.name()));
    }
}

auto CentralDifference::interpolate(
    double m_dot,
    const mesh::Cell& cell,
    const mesh::Cell& neighbor,
    const mesh::Face& face,
    const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme) -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    auto face_grad_phi = grad_scheme->gradient_at_face(face);
    auto d_Cf = face.center() - cell.center();
    auto a_C = std::max(m_dot, 0.0);
    auto b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    auto a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

auto Upwind::interpolate(
    double m_dot,
    const mesh::Cell& cell,                                           // NOLINT
    const mesh::Cell& neighbor,                                       // NOLINT
    const mesh::Face& face,                                           // NOLINT
    const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme) // NOLINT
    -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    auto a_C = std::max(m_dot, 0.0);
    // in case 'neighbor' is the upstream cell
    auto a_N = -std::max(-m_dot, 0.0);

    auto b = 0.0;

    return {a_C, a_N, b};
}

auto SecondOrderUpwind::interpolate(
    double m_dot,
    const mesh::Cell& cell,
    const mesh::Cell& neighbor,
    const mesh::Face& face,
    const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme) -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    const auto& face_grad_phi = grad_scheme->gradient_at_face(face);
    const auto& cell_grad_phi = grad_scheme->gradient_at_cell(cell);
    const auto& neighbor_grad_phi = grad_scheme->gradient_at_cell(neighbor);

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

auto QUICK::interpolate(double m_dot,
                        const mesh::Cell& cell,
                        const mesh::Cell& neighbor,
                        const mesh::Face& face,
                        const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
    -> CoeffsTriplet {
    // in case `cell` is the upstream cell
    const auto& face_grad_phi = grad_scheme->gradient_at_face(face);
    const auto& cell_grad_phi = grad_scheme->gradient_at_cell(cell);
    const auto& neighbor_grad_phi = grad_scheme->gradient_at_cell(neighbor);

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
