#include "cell.h"

namespace prism::mesh {
Cell::Cell(const std::vector<Face>& faces, std::vector<std::size_t>& faces_ids,
           std::size_t cell_id)
    : _id(cell_id) {
    Vector3d geo_center {0.0, 0.0, 0.0};

    // calculate geometric center
    for (auto face_id : faces_ids) {
        geo_center += faces[face_id].center();
    }

    geo_center /= faces_ids.size();

    // for each cell face construct a pyramid,
    // with the face as the pyramid base, and cell geometric center as the apex.
    for (auto face_id : faces_ids) {
        const auto& face = faces[face_id];
        const auto& face_center = face.center();
        const auto& face_area = face.area();

        Vector3d pyramid_center {(0.75 * face_center) + (0.25 * geo_center)};
        double pyramid_vol = ((1.0 / 3.0) * face_area) * (face_center - geo_center).norm();

        _volume += pyramid_vol;
        _center += (pyramid_vol * pyramid_center);
    }

    _center /= _volume;
    faces_ids = std::move(faces_ids);
}
} // namespace prism::mesh