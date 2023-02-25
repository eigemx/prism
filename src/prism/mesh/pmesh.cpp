#include "pmesh.h"

namespace prism::mesh {
// define all functions in PMesh class in pmesh.h
PMesh::PMesh(std::vector<Cell>&& cells, std::vector<Face>&& faces)
    : _cells(std::move(cells)), _faces(std::move(faces)) {}

PMesh::PMesh(std::vector<Cell>&& cells,
             std::vector<Face>&& faces,
             BoundaryConditions&& boundary_patches)
    : _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_conditions(std::move(boundary_patches)) {}

void PMesh::set_boundary_conditions(BoundaryConditions&& boundary_patches) {
    _boundary_conditions = std::move(boundary_patches);
}

void PMesh::set_boundary_conditions(const BoundaryConditions& boundary_patches) {
    _boundary_conditions = boundary_patches;
}

auto PMesh::non_ortho(std::size_t face_id) const -> double {
    return non_ortho(_faces.at(face_id));
}

auto PMesh::non_ortho(const Face& face) const -> double {
    if (!face.neighbor()) {
        throw std::runtime_error("PMesh::non_ortho() was called on a boundary face.");
    }

    const auto& owner_center = _cells[face.owner()].center();
    const auto& neigh_center = _cells[face.neighbor().value()].center();

    auto v = neigh_center - owner_center;
    auto v_norm = v.norm();

    return std::acos(v.dot(face.normal()) / v_norm) * (180. / pi);
}


} // namespace prism::mesh