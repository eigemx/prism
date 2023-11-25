#include "pmesh.h"

#include "../constants.h"

namespace prism::mesh {
PMesh::PMesh(std::vector<Vector3d> vertices,
             std::vector<Cell> cells,
             std::vector<Face> faces,
             std::vector<BoundaryPatch> boundary_patches,
             std::vector<std::size_t> boundary_faces_ids,
             std::vector<std::size_t> interior_faces_ids) noexcept
    : _vertices(std::move(vertices)),
      _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_patches(std::move(boundary_patches)),
      _boundary_faces_ids(std::move(boundary_faces_ids)),
      _interior_faces_ids(std::move(interior_faces_ids)),
      _n_cells(_cells.size()),
      _n_faces(_faces.size()) {
    // TODO: Check if inputs constitutes a valid mesh.
}

auto PMesh::face_non_ortho(std::size_t face_id) const -> double {
    return face_non_ortho(_faces.at(face_id));
}

auto PMesh::face_non_ortho(const Face& face) const -> double {
    if (!face.neighbor()) {
        throw std::runtime_error("PMesh::face_non_ortho() was called on a boundary face.");
    }

    const auto& owner_center = _cells[face.owner()].center();
    const auto& neigh_center = _cells[face.neighbor().value()].center();

    auto v = neigh_center - owner_center;
    auto v_norm = v.norm();

    return std::acos(v.dot(face.normal()) / (v_norm + EPSILON)) * (180. / PI);
}

// TODO: face_boundary_patch() methods don't check if face is boundary or not
// this is to avoid branching in the code, but it might be better to check
// think this over
auto PMesh::face_boundary_patch(std::size_t face_id) const -> const BoundaryPatch& {
    return face_boundary_patch(_faces[face_id]);
}

auto PMesh::face_boundary_patch(const Face& face) const -> const BoundaryPatch& {
    return _boundary_patches[face.boundary_patch_id().value()];
}


} // namespace prism::mesh