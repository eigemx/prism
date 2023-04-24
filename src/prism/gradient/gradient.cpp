#include "gradient.h"

#include "../print.h"

namespace prism::gradient {
auto boundary_face_phi(const mesh::Cell& cell, const mesh::Face& f, const ScalarField& field)
    -> std::optional<double> {
    /** @brief Calculates the value of the field at the boundary face
     *
     * @param cell Cell containing the face
     * @param f Boundary face
     * @param field Field to calculate the value
     * @return std::optional<double> Value of the field at the boundary face
     */

    const auto& mesh = field.mesh();

    auto face_boundary_patch_id = f.boundary_patch_id().value();
    const auto& boundary_patch = mesh.boundary_patches()[face_boundary_patch_id];

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty: {
            return std::nullopt;
        }

        case mesh::BoundaryPatchType::Fixed: {
            const auto& phi_name = field.name();
            return boundary_patch.get_scalar_bc(phi_name);
        }

        case mesh::BoundaryPatchType::Symmetry: {
            return field.data()[cell.id()];
        }

        default: {
            throw std::runtime_error(format("Non-Implemented boundary type for boundary patch {}",
                                            boundary_patch.name()));
        }
    }

    return std::nullopt;
}

auto GreenGauss::gradient(const mesh::Cell& cell, const ScalarField& field) -> Vector3d {
    /** @brief Calculates the gradient of a scalar field using Green Gauss theorem
     *
     * @param cell Cell to calculate the gradient at its center
     * @param field Field to calculate the gradient
     * @return Vector3d Gradient of the field at the cell center
     */

    Vector3d grad {0., 0., 0.};
    const auto& mesh = field.mesh();

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.faces()[face_id];
        auto Sf = face.area_vector();

        if (face.is_boundary()) {
            auto phi_opt = boundary_face_phi(cell, face, field);

            // skip empty boundary faces
            if (!phi_opt.has_value()) {
                continue;
            }

            auto face_phi = phi_opt.value();
            grad += Sf * face_phi;
            continue;
        }

        auto neighbor_cell_id =
            (face.owner() == cell.id()) ? face.neighbor().value() : face.owner();
        const auto& neighbor_cell = mesh.cells()[neighbor_cell_id];

        // gradient without skewness correction
        // TODO: implement skewness correction
        auto gc = mesh::PMesh::cells_weighting_factor(cell, neighbor_cell, face);
        auto face_phi = gc * field.data()[cell.id()];
        face_phi += (1 - gc) * field.data()[neighbor_cell_id];

        grad += Sf * face_phi;
    }

    return grad / cell.volume();
}

} // namespace prism::gradient