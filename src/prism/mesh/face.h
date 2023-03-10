#pragma once

#include <optional>
#include <tuple>
#include <vector>

#include "../types.h"

namespace prism::mesh {

class Face {
  public:
    Face(const std::vector<Vector3d>& face_vertices,
         std::vector<std::size_t> face_vertices_ids) noexcept;

    auto inline area() const noexcept -> const double& { return _area; }
    auto inline normal() const noexcept -> const Vector3d& { return _normal; }
    auto inline center() const noexcept -> const Vector3d& { return _center; }

    auto inline id() const noexcept -> std::size_t { return _id; }
    void inline set_id(std::size_t face_id) noexcept { _id = face_id; }
    auto inline vertices_ids() const noexcept -> const std::vector<std::size_t>& {
        return _vertices_ids;
    }

    auto inline owner() const -> std::size_t { return _owner.value(); }
    void inline set_owner(std::size_t owner_id) noexcept { _owner = owner_id; }

    // This function is added mainly for debugging purposes,
    // because ach face must have an owner, so this will return false if something is wrong
    auto inline has_owner() const -> bool { return _owner.has_value(); }

    auto inline neighbor() const -> const std::optional<std::size_t>& { return _neighbor; }
    auto inline has_neighbor() const -> bool { return _neighbor.has_value(); }
    auto inline set_neighbor(std::size_t nei_id) noexcept { _neighbor = nei_id; }

    auto inline boundary_patch_id() const -> const std::optional<std::size_t>& {
        return _boundary_patch_id;
    }
    void inline set_boundary_patch_id(std::size_t boundary_patch_id) noexcept {
        _boundary_patch_id = boundary_patch_id;
    }

  private:
    void set_face_attributes(const std::vector<Vector3d>& vertices);

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
