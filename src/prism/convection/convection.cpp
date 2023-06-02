#include "convection.h"

namespace prism::convection {

auto inline face_flow_rate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

void CentralDifference::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();

    // get adjacent cell sharing `face` with `cell`
    auto adjacent_cell_id = face.is_owned_by(cell_id) ? face.neighbor().value() : face.owner();

    const auto& adj_cell = _mesh.cell(adjacent_cell_id);

    // face area vector
    const auto& S_f = face.area_vector();

    // interpolated velocity vector at face centroid
    const auto& U_f = 0.5 * (_U[cell_id] + _U[adjacent_cell_id]);

    auto m_dot_f = face_flow_rate(_rho, U_f, S_f);

    coeff_matrix().coeffRef(cell_id, cell_id) += m_dot_f / 2;
    coeff_matrix().coeffRef(cell_id, adjacent_cell_id) += m_dot_f / 2;
}


void CentralDifference::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        case mesh::BoundaryPatchType::Fixed: {
        }

        default:
            throw std::runtime_error(
                format("convection::CentralDifference::apply_boundary(): "
                       "Non-implemented boundary type for boundary patch: '{}'",
                       boundary_patch.name()));
    }
}


} // namespace prism::convection
