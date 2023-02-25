#pragma once

#include <optional>
#include <tuple>
#include <vector>

#include "../types.h"

namespace prism::mesh {

class Face {
public:
    Face(std::vector<Vector3d>&& face_vertices) noexcept;
    Face(const std::vector<Vector3d>& face_vertices) noexcept;

    auto inline area() const noexcept -> const double& { return _area; }
    auto inline normal() const noexcept -> const Vector3d& { return _normal; }
    auto inline center() const noexcept -> const Vector3d& { return _center; }
    auto aspect_ratio() const noexcept -> double;

    auto inline id() const noexcept -> std::size_t { return _id; }
    void inline set_id(std::size_t face_id) noexcept { _id = face_id; }

    auto inline owner() const -> std::size_t { return _owner.value(); }
    void inline set_owner(std::size_t owner_id) { _owner = owner_id; }

    // This function is added mainly for debugging purposes,
    // because ach face must have an owner, so this will return false if something is wrong
    auto inline has_owner() const -> bool { return _owner.has_value(); }

    auto inline neighbor() const -> const std::optional<std::size_t>& { return _neighbor; }
    auto inline has_neighbor() const -> bool { return _neighbor.has_value(); }
    auto inline set_neighbor(std::size_t nei_id) noexcept { _neighbor = nei_id; }

private:
    void set_face_attributes();

    std::size_t vertices_count {0};
    std::vector<Vector3d> _vertices;
    double _area {0.0};
    Vector3d _normal {0.0, 0.0, 0.0};
    Vector3d _center {0.0, 0.0, 0.0};
    std::size_t _id {0};
    std::optional<std::size_t> _owner {std::nullopt};
    std::optional<std::size_t> _neighbor {std::nullopt};
};
} // namespace prism::mesh
