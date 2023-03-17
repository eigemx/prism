#include "gradient.h"

#include "../print.h"

namespace prism::gradient {
auto boundary_face_phi(const mesh::Cell& cell, const mesh::Face& f, const ScalarField& field)
    -> std::optional<double> {
    const auto& mesh = field.mesh();

    auto face_boundary_patch_id = f.boundary_patch_id().value();
    const auto& boundary_patch = mesh.boundary_patches()[face_boundary_patch_id];

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty: {
            return std::nullopt;
        }

        case mesh::BoundaryPatchType::Fixed: {
            const auto& phi_name = field.name();
            return boundary_patch.get_scalar_attribute(phi_name);
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
    auto grad = Vector3d::Zero();
    const auto& mesh = field.mesh();

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.faces()[face_id];
        auto Sf = face.normal() * face.area();

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

    grad /= cell.volume();
    return grad;
}

} // namespace prism::gradient