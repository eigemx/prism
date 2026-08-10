#pragma once

#include <optional>
#include <vector>

#include "prism/types.h"

namespace prism::mesh {

class Face {
  public:
    Face(const std::vector<Vector3d>& face_vertices,
         std::vector<size_t> face_vertices_ids) noexcept;

    auto area() const noexcept -> double { return _area; }
    auto normal() const noexcept -> const Vector3d& { return _normal; }
    auto center() const noexcept -> const Vector3d& { return _center; }
    auto areaVector() const noexcept -> Vector3d { return _area * _normal; }

    auto id() const noexcept -> size_t { return _id; }

    /// TODO: This should be setId() instead of id(), to make the purpose of the operation clear
    auto id() noexcept -> size_t& { return _id; }
    auto verticesIds() const noexcept -> const std::vector<size_t>& { return _vertices_ids; }

    auto owner() const -> size_t { return _owner.value(); }
    void setOwner(size_t owner_id) noexcept { _owner = owner_id; }

    // This function is added mainly for debugging purposes,
    // because each face must have an owner, so this will return false if something is wrong
    auto hasOwner() const -> bool { return _owner.has_value(); }
    auto isOwnedBy(const size_t& cell_id) const -> bool { return _owner.value() == cell_id; }

    auto neighbor() const -> const std::optional<size_t>& { return _neighbor; }
    auto isInterior() const -> bool { return _neighbor.has_value(); }
    auto setNeighbor(size_t nei_id) noexcept { _neighbor = nei_id; }
    auto isBoundary() const -> bool { return !isInterior(); }

    auto boundaryPatchId() const -> const std::optional<size_t>& { return _boundary_patch_id; }
    void setBoundaryPatchId(size_t boundary_patch_id) noexcept {
        _boundary_patch_id = boundary_patch_id;
    }

    /// True if this face belongs to an empty boundary patch. Interior faces and faces on
    /// non-empty patches always return false.
    auto isEmpty() const noexcept -> bool { return _is_empty; }
    void setIsEmpty(bool is_empty) noexcept { _is_empty = is_empty; }

  private:
    size_t _id {0};
    size_t _vertices_count {0};

    double _area {0.0};
    Vector3d _normal {0.0, 0.0, 0.0};
    Vector3d _center {0.0, 0.0, 0.0};

    std::vector<size_t> _vertices_ids;
    std::optional<size_t> _owner {std::nullopt};
    std::optional<size_t> _neighbor {std::nullopt};
    std::optional<size_t> _boundary_patch_id {std::nullopt};
    bool _is_empty {false};
};
} // namespace prism::mesh
