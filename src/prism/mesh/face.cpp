#include "face.h"


namespace prism::mesh {

void Face::set_face_attributes(const std::vector<Vector3d>& face_vertices, Vector3d& geo_center) {
    if (vertices_count == 3) {
        const auto& v0 = face_vertices[0];
        const auto& v1 = face_vertices[1];
        const auto& v2 = face_vertices[2];

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
        const Vector3d& v1 = face_vertices[i];
        const Vector3d& v2 = face_vertices[(i + 1) % vertices_count];
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

Face::Face(const std::vector<Vector3d>& face_vertices) : vertices_count(face_vertices.size()) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < vertices_count; ++i) {
        geo_center += face_vertices[i];
    }

    // Calculate face geometric center.
    geo_center /= vertices_count;

    set_face_attributes(face_vertices, geo_center);
}

} // namespace prism::mesh
