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

auto internal_face_phi(const mesh::Cell& cell, const mesh::Face& f) -> double {}

GreenGauss::GreenGauss(const ScalarField& field) : _field(field) {
    _cell_gradients.resize(field.mesh().cells().size());
}

auto GreenGauss::gradient(const mesh::Cell& cell) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = _field.mesh();

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.faces()[face_id];
        auto Sf = face.area_vector();

        if (face.is_boundary()) {
            const auto& face_phi = boundary_face_phi(cell, face, _field);

            // skip empty boundary faces
            if (!face_phi.has_value()) {
                continue;
            }

            grad += Sf * face_phi.value();
            continue;
        }

        auto neighbor_cell_id =
            (face.owner() == cell.id()) ? face.neighbor().value() : face.owner();
        const auto& neighbor_cell = mesh.cells()[neighbor_cell_id];

        // gradient without skewness correction
        // TODO: implement skewness correction
        auto gc = mesh::PMesh::cells_weighting_factor(cell, neighbor_cell, face);

        // phi_f = gc * phi_c + (1 - gc) * phi_n
        auto face_phi = gc * _field.data()[cell.id()];
        face_phi += (1 - gc) * _field.data()[neighbor_cell_id];

        grad += Sf * face_phi;
    }

    return grad / cell.volume();
}

} // namespace prism::gradient