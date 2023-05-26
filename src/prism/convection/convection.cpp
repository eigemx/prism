#include "convection.h"


namespace prism::convection {
auto inline mass_flow_rate(double rho, const Vector3d& U, const Vector3d& Sf) -> double {
    return rho * U.dot(Sf);
}

void Upwind::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    // get neighbor cell
    bool is_owner = face.owner() == cell.id();
    auto cell_id = cell.id();
    const auto& adj_cell_id = is_owner ? face.neighbor().value() : face.owner();
    const auto& adj_cell = _mesh.cells()[adj_cell_id];

    // calculate mass flow rate across face
    const auto& U_c = _U[cell_id];
    const auto& U_n = _U[adj_cell_id];

    // assign mass flow rate to the max value of the two cells velocities
    const auto& Sf = face.area_vector();
    auto m_dot_e = std::max(mass_flow_rate(_rho, U_c, Sf), mass_flow_rate(_rho, U_n, Sf));
}
} // namespace prism::convection