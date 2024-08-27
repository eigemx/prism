#pragma once

#include <optional>
#include <vector>

#include "../types.h"

namespace prism::mesh {

class Face {
  public:
    Face(const std::vector<Vector3d>& face_vertices,
         std::vector<std::size_t> face_vertices_ids) noexcept;

    auto inline area() const noexcept -> double { return _area; }
    auto inline normal() const noexcept -> const Vector3d& { return _normal; }
    auto inline center() const noexcept -> const Vector3d& { return _center; }
    auto inline areaVector() const noexcept -> Vector3d { return _area * _normal; }

    auto inline id() const noexcept -> std::size_t { return _id; }
    auto inline id() noexcept -> std::size_t& { return _id; }
    auto inline verticesIds() const noexcept -> const std::vector<std::size_t>& {
        return _vertices_ids;
    }

    auto inline owner() const -> std::size_t { return _owner.value(); }
    void inline setOwner(std::size_t owner_id) noexcept { _owner = owner_id; }

    // This function is added mainly for debugging purposes,
    // because each face must have an owner, so this will return false if something is wrong
    auto inline hasOwner() const -> bool { return _owner.has_value(); }
    auto inline isOwnedBy(const std::size_t& cell_id) const -> bool {
        return _owner.value() == cell_id;
    }

    auto inline neighbor() const -> const std::optional<std::size_t>& { return _neighbor; }
    auto inline isInterior() const -> bool { return _neighbor.has_value(); }
    auto inline setNeighbor(std::size_t nei_id) noexcept { _neighbor = nei_id; }
    auto inline isBoundary() const -> bool { return !isInterior(); }

    auto inline boundaryPatchId() const -> const std::optional<std::size_t>& {
        return _boundary_patch_id;
    }
    void inline setBoundaryPatchId(std::size_t boundary_patch_id) noexcept {
        _boundary_patch_id = boundary_patch_id;
    }

  private:
    std::size_t _id {0};
    std::size_t _vertices_count {0};

    double _area {0.0};
    Vector3d _normal {0.0, 0.0, 0.0};
    Vector3d _center {0.0, 0.0, 0.0};

    std::vector<std::size_t> _vertices_ids;
    std::optional<std::size_t> _owner {std::nullopt};
    std::optional<std::size_t> _neighbor {std::nullopt};
    std::optional<std::size_t> _boundary_patch_id {std::nullopt};
};
} // namespace prism::mesh
