#include "face.h"


namespace prism::mesh {

void Face::set_face_attributes(const std::vector<Vector3d>& face_vertices, Vector3d& geo_center) {
    // Calculate face geometric center.
    geo_center = geo_center / vertices_count;

    if (vertices_count == 3) {
        // calculate area and normal vector of triangle
        _normal =
            (face_vertices[1] - face_vertices[0]).cross(face_vertices[2] - face_vertices[0]);
        _area = _normal.norm() / 2.0;
        _center = (face_vertices[0] + face_vertices[1] + face_vertices[2]) / 3.0;
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

Face::Face(const std::vector<std::size_t>& face,
           const std::vector<std::array<double, 3>>& vertices)
    : vertices_count(face.size()) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    auto face_vertices = std::vector<Vector3d>();
    face_vertices.reserve(vertices_count);

    for (std::size_t i = 0; i < vertices_count; ++i) {
        const auto& vertex = vertices[face[i]];
        face_vertices.emplace_back(vertex[0], vertex[1], vertex[2]);
        geo_center += face_vertices[i];
    }

    set_face_attributes(face_vertices, geo_center);
}

Face::Face(const std::vector<std::size_t>& face, const std::vector<Vector3d>& vertices)
    : vertices_count(face.size()) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    auto face_vertices = std::vector<Vector3d>();
    face_vertices.reserve(vertices_count);

    for (std::size_t i = 0; i < vertices_count; ++i) {
        const auto& vertex = vertices[face[i]];
        face_vertices.emplace_back(vertex);
        geo_center += face_vertices[i];
    }

    set_face_attributes(face_vertices, geo_center);
}
} // namespace prism::mesh