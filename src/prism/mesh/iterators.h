#pragma once

#include <span>

#include "face.h"

namespace prism::mesh::iterators {
struct FaceIterator {
    using iterator_category = std::input_iterator_tag;
    using value_type = Face;
    using difference_type = std::ptrdiff_t;
    using pointer = Face*;
    using reference = Face&;

    FaceIterator(std::span<const Face> faces,
                 std::span<const std::size_t> ids,
                 std::size_t position);

    auto operator++() -> FaceIterator&;
    auto operator++(int) -> FaceIterator;
    auto operator*() const -> const Face&;
    auto operator->() const -> const Face*;
    auto operator==(const FaceIterator& other) const -> bool;
    auto operator!=(const FaceIterator& other) const -> bool;

  private:
    std::span<const Face> _faces;
    std::span<const std::size_t> _ids;
    std::size_t _current {};
};

struct BoundaryFaces {
    BoundaryFaces(std::span<const Face> faces,
                  std::span<const std::size_t> boundary_faces_ids);

    auto begin() const -> iterators::FaceIterator;
    auto end() const -> iterators::FaceIterator;

  private:
    std::span<const Face> _faces;
    std::span<const std::size_t> _boundary_faces_ids;
};

struct InteriorFaces {
    InteriorFaces(std::span<const Face> faces,
                  std::span<const std::size_t> interior_faces_ids);

    auto begin() const -> iterators::FaceIterator;
    auto end() const -> iterators::FaceIterator;

  private:
    std::span<const Face> _faces;
    std::span<const std::size_t> _interior_faces_ids;
};

} // namespace prism::mesh::iterators