#include <cmath>
#include <cstddef>

#include "gradient.h"
#include "prism/mesh/utilities.h"

namespace prism::gradient {

auto GreenGauss::correctSkewness(const mesh::Face& face,
                                 const mesh::Cell& cell,
                                 const mesh::Cell& nei) const -> double {
    auto grad_sum = _cell_gradients[cell.id()] + _cell_gradients[nei.id()];
    auto vec = face.center() - (0.5 * (cell.center() + nei.center()));

    return 0.5 * grad_sum.dot(vec);
}

GreenGauss::GreenGauss(field::Scalar field) : _field(field), IGradient(field) { // NOLINT
    // We need to perform a first run for calculating gradient at cells,
    // to make the cell gradient vector _cell_gradients available if the user desires to call
    // gradient_at_face() which requires a first run of gradient calculations, to perform
    // correction for faces with skewness.
    const std::size_t n_cells = _field.mesh().nCells();


    _cell_gradients.reserve(n_cells);
    for (const auto& cell : _field.mesh().cells()) {
        // caclulate the gradient without skewness correction
        _cell_gradients.emplace_back(gradAtCell_(cell, false));
    }
}

auto GreenGauss::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    return gradAtCell_(cell, true);
}

auto GreenGauss::gradAtCell_(const mesh::Cell& cell, bool correct_skewness) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = _field.mesh();

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.face(face_id);

        // This is a boundary face
        if (face.is_boundary()) {
            grad += boundaryFaceIntegral(face);
            continue;
        }

        // This is an internal face
        // Area normal vector, poitning out of the cell
        auto Sf = mesh::outward_area_vector(face, cell);
        const auto& nei = _field.mesh().otherSharingCell(cell, face);
        auto face_phi = 0.5 * (_field[cell.id()] + _field[nei.id()]);

        if (correct_skewness) {
            face_phi += correctSkewness(face, cell, nei);
        }
        grad += Sf * face_phi;
    }
    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients[cell.id()] = grad;

    return grad;
}

auto GreenGauss::boundaryFaceIntegral(const mesh::Face& face) -> Vector3d {
    // returns the Green-Gauss face integral at boundary face `face`
    const auto& boundary_patch = _field.mesh().boundary_patch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(_field.name());

    if (boundary_patch.isEmpty()) {
        return {0.0, 0.0, 0.0};
    }

    auto phi = _field.valueAtFace(face);
    return phi * face.area_vector();
}
} // namespace prism::gradient