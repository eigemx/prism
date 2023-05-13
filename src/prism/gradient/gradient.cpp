#include "gradient.h"

#include <cassert>
#include <iostream>

#include "../print.h"

namespace prism::gradient {


auto inline boundary_face_phi(const mesh::Cell& cell,
                              const mesh::Face& f,
                              const ScalarField& field) -> std::optional<double> {
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

auto inline skewness_correction(const mesh::Cell& C,
                                const mesh::Cell& F,
                                const mesh::Face& f,
                                const std::vector<Vector3d>& grad_field) -> double {
    auto grad_sum = grad_field[C.id()] + grad_field[F.id()];
    auto vec = f.center() - (0.5 * (C.center() + F.center()));

    return 0.5 * grad_sum.dot(vec);
}

GreenGauss::GreenGauss(const ScalarField& field) : _field(field) {
    _cell_gradients.resize(field.mesh().cells().size());

    for (auto& cell_gradient : _cell_gradients) {
        cell_gradient = Vector3d {0., 0., 0.};
    }
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

            // TODO: this only works for fixed boundary conditions
            grad += Sf * face_phi.value();
            continue;
        }

        // This is an internal face
        bool is_cell_owner = face.owner() == cell.id();
        auto neighbor_cell_id = is_cell_owner ? face.neighbor().value() : face.owner();

        // update normal vector to always be pointing out of the cell
        if (!is_cell_owner) {
            Sf *= -1;
        }

        const auto& neighbor_cell = mesh.cells()[neighbor_cell_id];

        auto gc = mesh::PMesh::cells_weighting_factor(cell, neighbor_cell, face);

        // phi_f = gc * phi_c + (1 - gc) * phi_n (Page 278)
        //auto face_phi = (gc * _field[cell.id()]) + ((1 - gc) * _field[neighbor_cell_id]);
        auto face_phi = 0.5 * (_field[cell.id()] + _field[neighbor_cell_id]);

        // skewness correction
        face_phi += skewness_correction(cell, neighbor_cell, face, _cell_gradients);

        grad += Sf * face_phi;
    }

    grad /= cell.volume();
    _cell_gradients[cell.id()] = grad;

    return grad;
}

} // namespace prism::gradient