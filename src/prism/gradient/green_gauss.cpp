#include <cassert>
#include <cmath>

#include "../print.h"
#include "gradient.h"

namespace prism::gradient {

auto inline skewness_correction(const mesh::Cell& C,
                                const mesh::Cell& F,
                                const mesh::Face& f,
                                const MatrixX3d& grad_field) -> double {
    auto grad_sum = grad_field.row(C.id()) + grad_field.row(F.id());
    auto vec = f.center() - (0.5 * (C.center() + F.center()));

    return 0.5 * grad_sum.dot(vec);
}

auto inline gradient_at_interior_face(const mesh::Cell& cell,
                                      const mesh::Face& face,
                                      const mesh::PMesh& mesh,
                                      const ScalarField& field,
                                      const MatrixX3d& grad_field) -> Vector3d {
    bool is_owned = face.is_owned_by(cell.id());
    auto neighbor_cell_id = is_owned ? face.neighbor().value() : face.owner();

    // update normal vector to always be pointing out of the cell
    auto Sf = face.area_vector() * std::pow(-1., static_cast<int>(is_owned));

    const auto& neighbor_cell = mesh.cell(neighbor_cell_id);

    auto face_phi = 0.5 * (field[cell.id()] + field[neighbor_cell_id]);

    // skewness correction
    face_phi += skewness_correction(cell, neighbor_cell, face, grad_field);

    return Sf * face_phi;
}

auto GreenGauss::gradient_at_cell(const mesh::Cell& cell) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = _field.mesh();

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.face(face_id);
        auto Sf = face.area_vector();

        // This is a boundary face
        if (face.is_boundary()) {
            grad += gradient_at_boundary_face(face, _field);
            continue;
        }

        // This is an internal face
        grad += gradient_at_interior_face(cell, face, mesh, _field, _cell_gradients);
    }

    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients.row(cell.id()) = grad;

    return grad;
}

auto GreenGauss::gradient_at_face(const mesh::Face& face) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = _field.mesh();

    if (face.is_boundary()) {
        return gradient_at_boundary_face(face, _field);
    }

    const auto& owner_cell = mesh.cell(face.owner());
    return gradient_at_interior_face(owner_cell, face, mesh, _field, _cell_gradients);
}

auto GreenGauss::gradient_field() -> VectorField {
    auto grad_field_name = _field.name() + "_grad";
    VectorField grad_field {grad_field_name, _field.mesh()};

    for (const auto& cell : _field.mesh().cells()) {
        grad_field.data().row(cell.id()) = gradient_at_cell(cell);
    }

    return grad_field;
}

} // namespace prism::gradient