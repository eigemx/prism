#include "face.h"


namespace prism::mesh {

void Face::set_face_attributes(const std::vector<Vector3d>& face_vertices, Vector3d& geo_center) {
    // Calculate face geometric center.
    geo_center = geo_center / vertices_count;

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
    _normal /= 2.0;
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