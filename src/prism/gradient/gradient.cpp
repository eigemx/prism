#include "gradient.h"

namespace prism::gradient {
auto GreenGauss::gradient(const mesh::Cell& c, const ScalarField& field) -> Vector3d {
    auto grad = Vector3d::Zero();
    const auto& mesh = field.mesh();

    for (auto face_id : c.faces_ids()) {
        std::size_t face_phi {};
        const auto& face = mesh.faces()[face_id];

        if (!face.has_neighbor()) {
            face_phi = face.owner();
            continue;
        }

        // get the other cell sharing the face
        std::size_t adj_cell_id {};

        if (face.owner() == c.id()) {
            //adj_cell_id = face.neighbor();
        } else {
            adj_cell_id = face.owner();
        }
    }
}
} // namespace prism::gradient