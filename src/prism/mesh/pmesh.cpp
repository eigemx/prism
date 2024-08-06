#include "pmesh.h"

#include <cassert>

#include "spdlog/spdlog.h"

namespace prism::mesh {

namespace detail {
FaceIterator::FaceIterator(const std::vector<Face>& faces,
                           const std::vector<std::size_t>& ids,
                           std::size_t position)
    : _faces(faces), _ids(ids), _current(position) {}

auto FaceIterator::operator++() -> FaceIterator& {
    _current++;
    return *this;
}

auto FaceIterator::operator++(int) -> FaceIterator {
    auto tmp = *this;
    ++(*this);
    return tmp;
}

auto FaceIterator::operator*() const -> const Face& {
    return _faces[_ids[_current]];
}
auto FaceIterator::operator->() const -> const Face* {
    return &_faces[_ids[_current]];
}

auto FaceIterator::operator==(const FaceIterator& other) const -> bool {
    return _current == other._current;
}

auto FaceIterator::operator!=(const FaceIterator& other) const -> bool {
    return !(*this == other);
}

BoundaryFaces::BoundaryFaces(const std::vector<Face>& faces,
                             const std::vector<std::size_t>& boundary_faces_ids)
    : _faces(faces), _boundary_faces_ids(boundary_faces_ids) {}

auto BoundaryFaces::begin() const -> detail::FaceIterator {
    return {_faces, _boundary_faces_ids, 0};
}

auto BoundaryFaces::end() const -> detail::FaceIterator {
    return {_faces, _boundary_faces_ids, _boundary_faces_ids.size()};
}

InteriorFaces::InteriorFaces(const std::vector<Face>& faces,
                             const std::vector<std::size_t>& interior_faces_ids)
    : _faces(faces), _interior_faces_ids(interior_faces_ids) {}

auto InteriorFaces::begin() const -> detail::FaceIterator {
    return {_faces, _interior_faces_ids, 0};
}
auto InteriorFaces::end() const -> detail::FaceIterator {
    return {_faces, _interior_faces_ids, _interior_faces_ids.size()};
}
} // namespace detail

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
    spdlog::debug("PMesh object created with {} cells, {} faces and {} vertices.",
                  _n_cells,
                  _n_faces,
                  _vertices.size());
    spdlog::debug("PMesh object has {} internal faces and {} boundary faces. ",
                  _interior_faces_ids.size(),
                  _boundary_faces_ids.size());

    _cells_volume.resize(_n_cells);
    for (const auto& cell : _cells) {
        _cells_volume[cell.id()] = cell.volume();
    }
}

// TODO: face_boundary_patch() methods don't check if face is boundary or not this is to avoid
// branching in the code, but it might be better to check think this over
auto PMesh::faceBoundaryPatch(std::size_t face_id) const -> const BoundaryPatch& {
    assert(
        face_id < _faces.size() &&
        "PMesh::face_boundary_patch() was called on a face with an index larger than mesh faces "
        "count");
    return faceBoundaryPatch(_faces[face_id]);
}

auto PMesh::faceBoundaryPatch(const Face& face) const -> const BoundaryPatch& {
    assert(face.is_boundary() && "PMesh::face_boundary_patch() was called on an interior face");
    return _boundary_patches[face.boundary_patch_id().value()];
}

auto PMesh::otherSharingCell(const Cell& c, const Face& f) const -> const Cell& {
    assert(f.is_interior() && "PMesh::other_sharing_cell() called on a boundary face!");
    auto n_id = f.owner() == c.id() ? f.neighbor().value() : f.owner();
    return _cells[n_id];
}


} // namespace prism::mesh