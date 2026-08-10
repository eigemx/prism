#include "pmesh.h"

#include <cassert>
#include <span>
#include <stdexcept>

#include "prism/log.h"

namespace prism::mesh {
PMesh::PMesh(std::vector<Vector3d> vertices,
             std::vector<Cell> cells,
             std::vector<Face> faces,
             std::vector<BoundaryPatch> boundary_patches,
             std::vector<FieldInfo> field_infos,
             std::vector<size_t> boundary_faces_ids,
             std::vector<size_t> interior_faces_ids) noexcept
    : _vertices(std::move(vertices)),
      _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_patches(std::move(boundary_patches)),
      _field_infos(std::move(field_infos)),
      _n_cells(_cells.size()),
      _n_faces(_faces.size()) {
    auto perm = buildFacePermutation(boundary_faces_ids, interior_faces_ids);

    reorderFaces(perm.old_to_new);
    remapFaceIds(perm.old_to_new);

    _n_boundary_faces = perm.n_boundary;
    _n_nonempty_boundary_faces = perm.n_nonempty;

    // Make each boundary face aware of whether it belongs to an empty patch, so callers can
    // skip empty faces with a single flag load instead of probing the owning patch.
    for (size_t face_id = 0; face_id < _n_boundary_faces; ++face_id) {
        const auto& face = _faces[face_id];
        _faces[face_id].setIsEmpty(_boundary_patches[face.boundaryPatchId().value()].isEmpty());
    }

    std::vector<size_t> {}.swap(perm.old_to_new);
    std::vector<size_t> {}.swap(boundary_faces_ids);
    std::vector<size_t> {}.swap(interior_faces_ids);

    /// TODO: Check if inputs constitutes a valid mesh.
    log::debug("prism::mesh::PMesh() object created with {} cells, {} faces and {} vertices.",
               _n_cells,
               _n_faces,
               _vertices.size());
    log::debug("prism::mesh::PMesh() object has {} internal faces and {} boundary faces. ",
               _n_faces - _n_boundary_faces,
               _n_boundary_faces);

    _cells_volume.resize(_n_cells);
    for (const auto& cell : _cells) {
        _cells_volume[cell.id()] = cell.volume();
    }
}

auto PMesh::buildFacePermutation(const std::vector<size_t>& boundary_faces_ids,
                                 const std::vector<size_t>& interior_faces_ids) const
    -> FacePermutation {
    auto n_faces = _faces.size();
    auto n_boundary = boundary_faces_ids.size();

    size_t n_nonempty = 0;
    for (auto old_id : boundary_faces_ids) {
        const auto& face = _faces[old_id];
        if (!_boundary_patches[face.boundaryPatchId().value()].isEmpty()) {
            ++n_nonempty;
        }
    }

    std::vector<size_t> old_to_new(n_faces);
    size_t nonempty_pos = 0;
    size_t empty_pos = n_nonempty;
    size_t interior_pos = n_boundary;

    for (auto old_id : boundary_faces_ids) {
        const auto& face = _faces[old_id];
        if (!_boundary_patches[face.boundaryPatchId().value()].isEmpty()) {
            old_to_new[old_id] = nonempty_pos++;
        } else {
            old_to_new[old_id] = empty_pos++;
        }
    }

    for (auto old_id : interior_faces_ids) {
        old_to_new[old_id] = interior_pos++;
    }

    return {
        .old_to_new = std::move(old_to_new), .n_boundary = n_boundary, .n_nonempty = n_nonempty};
}

void PMesh::reorderFaces(const std::vector<size_t>& old_to_new) {
    auto n_faces = _faces.size();

    std::vector<size_t> new_to_old(n_faces);
    for (size_t old_id = 0; old_id < n_faces; ++old_id) {
        new_to_old[old_to_new[old_id]] = old_id;
    }

    std::vector<Face> new_faces;
    new_faces.reserve(n_faces);
    for (size_t new_id = 0; new_id < n_faces; ++new_id) {
        auto old_id = new_to_old[new_id];
        new_faces.emplace_back(std::move(_faces[old_id]));
        new_faces.back().id() = new_id;
    }
    _faces = std::move(new_faces);
}

void PMesh::remapFaceIds(const std::vector<size_t>& old_to_new) {
    for (auto& cell : _cells) {
        auto& ids = cell.facesIds();
        for (auto& face_id : ids) {
            face_id = old_to_new[face_id];
        }
    }

    for (auto& patch : _boundary_patches) {
        auto& ids = patch.facesIds();
        for (auto& face_id : ids) {
            face_id = old_to_new[face_id];
        }
    }
}

auto PMesh::vertices() const noexcept -> const std::vector<Vector3d>& {
    return _vertices;
}

auto PMesh::cells() const noexcept -> const std::vector<Cell>& {
    return _cells;
}

auto PMesh::cells() noexcept -> std::vector<Cell>& {
    return _cells;
}

auto PMesh::cell(size_t cell_id) const -> const Cell& {
    return _cells[cell_id];
}

auto PMesh::cell(size_t cell_id) noexcept -> Cell& {
    return _cells[cell_id];
}

auto PMesh::faces() const noexcept -> const std::vector<Face>& {
    return _faces;
}
auto PMesh::faces() noexcept -> std::vector<Face>& {
    return _faces;
}

auto PMesh::face(size_t face_id) const -> const Face& {
    return _faces[face_id];
}

auto PMesh::face(size_t face_id) noexcept -> Face& {
    return _faces[face_id];
}

auto PMesh::boundaryPatches() const noexcept -> const std::vector<BoundaryPatch>& {
    return _boundary_patches;
}

auto PMesh::boundaryPatch(const Face& face) const noexcept -> const BoundaryPatch& {
    assert(face.isBoundary() && face.boundaryPatchId().has_value());
    return _boundary_patches[face.boundaryPatchId().value()];
}

auto PMesh::boundaryPatch(const std::string& name) const -> const BoundaryPatch& {
    for (const auto& patch : _boundary_patches) {
        if (patch.name() == name) {
            return patch;
        }
    }
    throw std::runtime_error(
        fmt::format("prism::mesh::PMesh::boundaryPatch(): Couldn't find boundary patch with name "
                    "`{}`",
                    name));
}

auto PMesh::faceBoundaryPatch(size_t face_id) const -> const BoundaryPatch& {
    assert(face_id < _faces.size() &&
           "prism::mesh::PMesh::faceBoundaryPatch() was called on a face with an index larger "
           "than mesh faces "
           "count");
    return faceBoundaryPatch(_faces[face_id]);
}

auto PMesh::faceBoundaryPatch(const Face& face) const -> const BoundaryPatch& {
    assert(face.isBoundary() &&
           "prism::mesh::PMesh::faceBoundaryPatch() was called on an interior face");
    return _boundary_patches[face.boundaryPatchId().value()];
}

auto PMesh::cellCount() const noexcept -> size_t {
    return _n_cells;
}

auto PMesh::faceCount() const noexcept -> size_t {
    return _n_faces;
}

auto PMesh::boundaryFaceCount() const noexcept -> size_t {
    return _n_boundary_faces;
}

auto PMesh::cellsVolumeVector() const noexcept -> const VectorXd& {
    return _cells_volume;
}

auto PMesh::otherSharingCell(const Cell& c, const Face& f) const -> const Cell& {
    /// TODO: This is dangerous if called for a boundary face (no neighbor). Fix.
    assert(f.isInterior() && "prism::mesh::Mesh::otherSharingCell() called on a boundary face!");
    auto n_id = f.owner() == c.id() ? f.neighbor().value() : f.owner();
    return _cells[n_id];
}

auto PMesh::interiorFaces() const noexcept -> std::span<const Face> {
    return {_faces.data() + _n_boundary_faces, _faces.size() - _n_boundary_faces};
}

auto PMesh::boundaryFaces() const noexcept -> std::span<const Face> {
    return {_faces.data(), _n_boundary_faces};
}

auto PMesh::nonEmptyBoundaryFaces() const noexcept -> std::span<const Face> {
    return {_faces.data(), _n_nonempty_boundary_faces};
}

auto PMesh::fieldsInfo() const noexcept -> const std::vector<FieldInfo>& {
    return _field_infos;
}


} // namespace prism::mesh
