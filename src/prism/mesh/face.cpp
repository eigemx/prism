#include "face.h"


namespace prism::mesh {

Face::Face(std::vector<Vector3d>&& face_vertices) noexcept
    : vertices_count(face_vertices.size()), _vertices(std::move(face_vertices)) {
    set_face_attributes();
}

Face::Face(const std::vector<Vector3d>& face_vertices) noexcept
    : vertices_count(face_vertices.size()), _vertices(face_vertices) {
    set_face_attributes();
}

void inline Face::set_face_attributes() {
    Vector3d geo_center {0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < vertices_count; ++i) {
        geo_center += _vertices[i];
    }

    // Calculate face geometric center.
    geo_center /= vertices_count;
    if (vertices_count == 3) {
        const auto& v0 = _vertices[0];
        const auto& v1 = _vertices[1];
        const auto& v2 = _vertices[2];

        // calculate area and normal vector of triangle
        _normal = (v1 - v0).cross(v2 - v0);
        _area = _normal.norm() / 2.0;
        _center = (v0 + v1 + v2) / 3.0;
        _normal /= _normal.norm();
        return;
    }

    // form triangular subfaces
    // each subface is constructed using an edge (from main face) as the base.
    // and main face geomertic center as its apex.
    for (std::size_t i = 0; i < vertices_count; ++i) {
        // set subface points
        const Vector3d& v1 = _vertices[i];
        const Vector3d& v2 = _vertices[(i + 1) % vertices_count];
        const Vector3d& v3 = geo_center;

        // calculate subface geometric center
        Vector3d subface_geo_center {(v1 + v2 + v3) / 3.0};

        // calculate area and normal vector of subface
        auto subface_normal {(v2 - v1).cross(v3 - v1)};
        auto subface_area {subface_normal.norm() / 2.0};

        _normal += subface_normal;
        _area += subface_area;

        _center += (subface_geo_center * subface_area);
    }

    _center /= _area;
    _normal /= _normal.norm();
}

auto Face::aspect_ratio() const noexcept -> double {
    // find the bounding box of the face, and calculate its aspect ratio.
    Vector3d min {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max()};

    Vector3d max {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                  std::numeric_limits<double>::lowest()};

    for (const auto& vertex : _vertices) {
        min = min.cwiseMin(vertex);
        max = max.cwiseMax(vertex);
    }

    Vector3d diff {max - min};
    return diff.maxCoeff() / (diff.minCoeff() + 1e-8);
}

} // namespace prism::mesh
