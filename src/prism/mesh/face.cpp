#include "face.h"

namespace prism::mesh {

Face::Face(const std::vector<Vector3d>& face_vertices,
           std::vector<std::size_t> vertices_ids) noexcept
    : _vertices_count(face_vertices.size()), _vertices_ids(std::move(vertices_ids)) {
    set_face_attributes(face_vertices);
}

void Face::set_face_attributes(const std::vector<Vector3d>& vertices) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < _vertices_count; ++i) {
        geo_center += vertices[i];
    }

    // Calculate face geometric center.
    geo_center /= _vertices_count;

    // form triangular subfaces
    // each subface is constructed using an edge (from main face) as the base.
    // and main face geomertic center as its apex.
    for (std::size_t i = 0; i < _vertices_count; ++i) {
        // set subface points
        const Vector3d& v1 = vertices[i];
        const Vector3d& v2 = vertices[(i + 1) % _vertices_count];
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

} // namespace prism::mesh
