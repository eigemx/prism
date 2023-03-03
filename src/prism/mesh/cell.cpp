#include "cell.h"

#include <cmath>
#include <iostream>

namespace prism::mesh {

Cell::Cell(const std::vector<Face>& faces,
           std::vector<std::size_t> faces_ids,
           std::vector<std::size_t> vertices_ids,
           std::size_t cell_id)
    : _id(cell_id), _faces_ids(std::move(faces_ids)), _vertices_ids(std::move(vertices_ids)) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    // calculate geometric center
    for (auto face_id : _faces_ids) {
        geo_center += faces[face_id].center();
    }

    geo_center /= _faces_ids.size();

    // for each cell face construct a pyramid,
    // with the face as the pyramid base, and cell geometric center as the apex.
    for (auto face_id : _faces_ids) {
        const auto& face = faces[face_id];
        const auto& face_center = face.center();
        const auto& face_area = face.area();
        const auto& face_normal = face.normal();

        double pyramid_vol {0.0};

        // calculate pyramid volume
        if (face.owner() == _id) {
            pyramid_vol = face_normal.dot(face_center - geo_center) * face_area;
        } else {
            pyramid_vol = face_normal.dot(geo_center - face_center) * face_area;
        }

        Vector3d pyramid_center {(face_center * 0.75) + (geo_center * 0.25)};

        _volume += pyramid_vol;
        _center += (pyramid_vol * pyramid_center);
    }

    _volume /= 3.0;
    _center /= _volume;
}

} // namespace prism::mesh