#include "pmesh.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace prism::mesh {

namespace iterators {
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

auto BoundaryFaces::begin() const -> iterators::FaceIterator {
    return {_faces, _boundary_faces_ids, 0};
}

auto BoundaryFaces::end() const -> iterators::FaceIterator {
    return {_faces, _boundary_faces_ids, _boundary_faces_ids.size()};
}

InteriorFaces::InteriorFaces(const std::vector<Face>& faces,
                             const std::vector<std::size_t>& interior_faces_ids)
    : _faces(faces), _interior_faces_ids(interior_faces_ids) {}

auto InteriorFaces::begin() const -> iterators::FaceIterator {
    return {_faces, _interior_faces_ids, 0};
}

auto InteriorFaces::end() const -> iterators::FaceIterator {
    return {_faces, _interior_faces_ids, _interior_faces_ids.size()};
}

} // namespace iterators

PMesh::PMesh(std::vector<Vector3d> vertices,
             std::vector<Cell> cells,
             std::vector<Face> faces,
             std::vector<BoundaryPatch> boundary_patches,
             std::vector<FieldInfo> field_infos,
             std::vector<std::size_t> boundary_faces_ids,
             std::vector<std::size_t> interior_faces_ids) noexcept
    : _vertices(std::move(vertices)),
      _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_patches(std::move(boundary_patches)),
      _field_infos(std::move(field_infos)),
      _boundary_faces_ids(std::move(boundary_faces_ids)),
      _interior_faces_ids(std::move(interior_faces_ids)),
      _n_cells(_cells.size()),
      _n_faces(_faces.size()) {
    // TODO: Check if inputs constitutes a valid mesh.
    spdlog::debug("prism::mesh::PMesh() object created with {} cells, {} faces and {} vertices.",
                  _n_cells,
                  _n_faces,
                  _vertices.size());
    spdlog::debug("prism::mesh::PMesh() object has {} internal faces and {} boundary faces. ",
                  _interior_faces_ids.size(),
                  _boundary_faces_ids.size());

    _cells_volume.resize(_n_cells);
    for (const auto& cell : _cells) {
        _cells_volume[cell.id()] = cell.volume();
    }

    std::copy_if(_boundary_faces_ids.begin(),
                 _boundary_faces_ids.end(),
                 std::back_inserter(_nonempty_boundary_faces_ids),
                 [this](const std::size_t& face_id) {
                     const auto& patch = faceBoundaryPatch(face_id);
                     return !patch.isEmpty();
                 });
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

auto PMesh::cell(std::size_t cell_id) const -> const Cell& {
    return _cells[cell_id];
}

auto PMesh::cell(std::size_t cell_id) noexcept -> Cell& {
    return _cells[cell_id];
}

auto PMesh::faces() const noexcept -> const std::vector<Face>& {
    return _faces;
}
auto PMesh::faces() noexcept -> std::vector<Face>& {
    return _faces;
}

auto PMesh::face(std::size_t face_id) const -> const Face& {
    return _faces[face_id];
}

auto PMesh::face(std::size_t face_id) noexcept -> Face& {
    return _faces[face_id];
}

auto PMesh::boundaryPatches() const noexcept -> const std::vector<BoundaryPatch>& {
    return _boundary_patches;
}
auto PMesh::boundaryPatch(const Face& face) const noexcept -> const BoundaryPatch& {
    assert(face.isBoundary() && face.boundaryPatchId().has_value());
    return _boundary_patches[face.boundaryPatchId().value()];
}

// TODO: faceBoundaryPatch() methods don't check if face is boundary or not this is to avoid
// branching in the code, but it might be better to check think this over
auto PMesh::faceBoundaryPatch(std::size_t face_id) const -> const BoundaryPatch& {
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

auto PMesh::cellCount() const noexcept -> std::size_t {
    return _n_cells;
}

auto PMesh::faceCount() const noexcept -> std::size_t {
    return _n_faces;
}

auto PMesh::boundaryFaceCount() const noexcept -> std::size_t {
    return _boundary_faces_ids.size();
}

auto PMesh::nonEmptyboundaryFaceCount() const noexcept -> std::size_t {
    return _nonempty_boundary_faces_ids.size();
}

auto PMesh::cellsVolumeVector() const noexcept -> const VectorXd& {
    return _cells_volume;
}

auto PMesh::otherSharingCell(const Cell& c, const Face& f) const -> const Cell& {
    assert(f.isInterior() && "prism::mesh::Mesh::otherSharingCell() called on a boundary face!");
    auto n_id = f.owner() == c.id() ? f.neighbor().value() : f.owner();
    return _cells[n_id];
}

auto PMesh::interiorFaces() const -> iterators::InteriorFaces {
    return {_faces, _interior_faces_ids};
}

// TODO: for boundaryFaces() and nonEmptyBoundaryFaces() we could iterate over the boundary
// patches instead, this allows us to get rid of _boundary_faces_ids and
// _nonempty_boundary_faces_ids vectors.
auto PMesh::boundaryFaces() const -> iterators::BoundaryFaces {
    return {_faces, _boundary_faces_ids};
}

auto PMesh::nonEmptyBoundaryFaces() const -> iterators::BoundaryFaces {
    return {_faces, _nonempty_boundary_faces_ids};
}

auto PMesh::fieldsInfo() const noexcept -> const std::vector<FieldInfo>& {
    return _field_infos;
}

PMeshPtr::PMeshPtr(const PMesh* ptr) : _ptr(ptr) {
    if (ptr == nullptr) {
        throw std::runtime_error("prism::mesh::PMeshPtr() was given a null pointer");
    }
}


} // namespace prism::mesh