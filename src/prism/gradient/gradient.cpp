#include "gradient.h"

#include <cassert>
#include <iostream>

#include "../print.h"

namespace prism::gradient {

auto inline skewness_correction(const mesh::Cell& C,
                                const mesh::Cell& F,
                                const mesh::Face& f,
                                const std::vector<Vector3d>& grad_field) -> double {
    auto grad_sum = grad_field[C.id()] + grad_field[F.id()];
    auto vec = f.center() - (0.5 * (C.center() + F.center()));

    return 0.5 * grad_sum.dot(vec);
}

auto inline boundary_face_gradient(const mesh::Face& face, const ScalarField& field) -> Vector3d {
    const auto& mesh = field.mesh();

    auto face_boundary_patch_id = face.boundary_patch_id().value();
    const auto& boundary_patch = mesh.boundary_patches()[face_boundary_patch_id];

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty:
        case mesh::BoundaryPatchType::Symmetry: {
            return Vector3d {0., 0., 0.};
        }

        case mesh::BoundaryPatchType::Fixed: {
            const auto& phi_name = field.name();
            auto phi = boundary_patch.get_scalar_bc(phi_name);
            return phi * face.area_vector();
        }

        case mesh::BoundaryPatchType::FixedGradient: {
            const auto& phi_name = field.name();
            auto flux = boundary_patch.get_scalar_bc(phi_name + "-flux");
            return flux * face.area_vector();
        }

        default:
            throw std::runtime_error(
                format("gradient/gradient.cpp boundary_face_gradient(): "
                       "Non-implemented boundary type for boundary patch: '{}'",
                       boundary_patch.name()));
    }
}

auto inline interior_face_gradient(const mesh::Cell& cell,
                                   const mesh::Face& face,
                                   const mesh::PMesh& mesh,
                                   const ScalarField& field,
                                   const std::vector<Vector3d>& grad_field) -> Vector3d {
    auto Sf = face.area_vector();
    bool is_cell_owner = face.owner() == cell.id();
    auto neighbor_cell_id = is_cell_owner ? face.neighbor().value() : face.owner();

    // update normal vector to always be pointing out of the cell
    if (!is_cell_owner) {
        Sf *= -1;
    }

    const auto& neighbor_cell = mesh.cells()[neighbor_cell_id];

    auto face_phi = 0.5 * (field[cell.id()] + field[neighbor_cell_id]);

    // skewness correction
    face_phi += skewness_correction(cell, neighbor_cell, face, grad_field);

    return Sf * face_phi;
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

        // This is a boundary face
        if (face.is_boundary()) {
            grad += boundary_face_gradient(face, _field);
            continue;
        }

        // This is an internal face
        grad += interior_face_gradient(cell, face, mesh, _field, _cell_gradients);
    }

    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients[cell.id()] = grad;

    return grad;
}

} // namespace prism::gradient