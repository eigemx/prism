#pragma once

#include <optional>
#include <tuple>
#include <vector>

#include "../types.h"

namespace prism::mesh {

class Face {
public:
    Face() = delete;
    Face(const std::vector<std::size_t>& face,
         const std::vector<std::array<double, 3>>& vertices);
    Face(const std::vector<std::size_t>& face, const std::vector<Vector3d>& vertices);

    auto inline area() const -> const double& { return _area; }
    auto inline normal() const -> const Vector3d& { return _normal; }
    auto inline center() const -> const Vector3d& { return _center; }

    auto inline id() const -> std::size_t { return _id; }
    void inline set_id(std::size_t face_id) { _id = face_id; }

    auto inline owner() const -> std::size_t { return _owner; }
    void inline set_owner(std::size_t owner_id) { _owner = owner_id; }

    auto inline neighbor() const -> const std::optional<std::size_t>& { return _neighbor; }
    auto inline set_neighbor(std::size_t nei_id) { _neighbor = nei_id; }

private:
    void inline set_face_attributes(const std::vector<Vector3d>& face_vertices,
                                    Vector3d& geo_center);
    std::size_t vertices_count;
    double _area {0.0};
    Vector3d _normal {0.0, 0.0, 0.0};
    Vector3d _center {0.0, 0.0, 0.0};
    std::size_t _id {0};
    std::size_t _owner {0};
    std::optional<std::size_t> _neighbor {std::nullopt};
};
} // namespace prism::mesh
