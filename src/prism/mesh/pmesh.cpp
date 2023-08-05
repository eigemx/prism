#include "pmesh.h"

#include "../constants.h"

namespace prism::mesh {
// define all functions in PMesh class in pmesh.h
PMesh::PMesh(std::vector<Vector3d> vertices,
             std::vector<Cell> cells,
             std::vector<Face> faces) noexcept
    : _vertices(std::move(vertices)),
      _cells(std::move(cells)),
      _faces(std::move(faces)),
      _n_cells(_cells.size()) {}

PMesh::PMesh(std::vector<Vector3d> vertices,
             std::vector<Cell> cells,
             std::vector<Face> faces,
             std::vector<BoundaryPatch> boundary_patches) noexcept
    : _vertices(std::move(vertices)),
      _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_patches(std::move(boundary_patches)),
      _n_cells(_cells.size()) {}

void PMesh::set_boundary_conditions(std::vector<BoundaryPatch> boundary_patches) {
    _boundary_patches = std::move(boundary_patches);
}

void PMesh::set_boundary_conditions(const std::vector<BoundaryPatch>& boundary_patches) {
    _boundary_patches = boundary_patches;
}

auto PMesh::face_non_ortho(std::size_t face_id) const -> double {
    return face_non_ortho(_faces.at(face_id));
}

auto PMesh::face_non_ortho(const Face& face) const -> double {
    if (!face.neighbor()) {
        throw std::runtime_error("PMesh::non_ortho() was called on a boundary face.");
    }

    const auto& owner_center = _cells[face.owner()].center();
    const auto& neigh_center = _cells[face.neighbor().value()].center();

    auto v = neigh_center - owner_center;
    auto v_norm = v.norm();

    return std::acos(v.dot(face.normal()) / (v_norm + EPSILON)) * (180. / PI);
}

auto PMesh::face_boundary_patch(std::size_t face_id) const -> const BoundaryPatch& {
    return face_boundary_patch(_faces.at(face_id));
}

auto PMesh::face_boundary_patch(const Face& face) const -> const BoundaryPatch& {
    if (face.neighbor()) {
        throw std::runtime_error("PMesh::boundary_patch() was called on an interior face.");
    }

    return _boundary_patches.at(face.boundary_patch_id().value());
}


} // namespace prism::mesh