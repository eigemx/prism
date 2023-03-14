#include "gradient.h"

namespace prism::gradient {
inline auto boundary_face_phi(const mesh::Face& face,
                              const VectorXd& phi_scalar_field,
                              const mesh::PMesh& mesh) -> double {
    const auto& boundary_patch = mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Wall:
            auto phi_wall =
                std::get<mesh::WallBoundaryData>(boundary_patch.data()).temperature.value();
    }
}
auto GreenGauss::gradient(const mesh::Cell& c,
                          const VectorXd& phi_scalar_field,
                          const mesh::PMesh& mesh,
                          const std::string& scalar_name) -> Vector3d {
    auto grad = Vector3d::Zero();

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