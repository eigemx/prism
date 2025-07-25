#include "iterators.h"

namespace prism::mesh::iterators {
FaceIterator::FaceIterator(std::span<const Face> faces,
                           std::span<const std::size_t> ids,
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

BoundaryFaces::BoundaryFaces(std::span<const Face> faces,
                             std::span<const std::size_t> boundary_faces_ids)
    : _faces(faces), _boundary_faces_ids(boundary_faces_ids) {}

auto BoundaryFaces::begin() const -> iterators::FaceIterator {
    return {_faces, _boundary_faces_ids, 0};
}

auto BoundaryFaces::end() const -> iterators::FaceIterator {
    return {_faces, _boundary_faces_ids, _boundary_faces_ids.size()};
}

InteriorFaces::InteriorFaces(std::span<const Face> faces,
                             std::span<const std::size_t> interior_faces_ids)
    : _faces(faces), _interior_faces_ids(interior_faces_ids) {}

auto InteriorFaces::begin() const -> iterators::FaceIterator {
    return {_faces, _interior_faces_ids, 0};
}

auto InteriorFaces::end() const -> iterators::FaceIterator {
    return {_faces, _interior_faces_ids, _interior_faces_ids.size()};
}

}