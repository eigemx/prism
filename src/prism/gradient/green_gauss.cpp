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


GreenGauss::GreenGauss(field::IScalar* field) : IGradient(field) { // NOLINT
    const std::size_t n_cells = field->mesh().nCells();

    _cell_gradients.reserve(n_cells);

    // TODO: replace this with std::transform
    for (const auto& cell : field->mesh().cells()) {
        _cell_gradients.emplace_back(Vector3d::Zero());
    }
}

auto GreenGauss ::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    return gradAtCell_(cell, true);
}

auto GreenGauss ::gradAtCellStored(const mesh::Cell& cell) -> Vector3d {
    return _cell_gradients[cell.id()];
}

auto GreenGauss ::gradAtCell_(const mesh::Cell& cell, bool correct_skewness) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = this->field()->mesh();

    for (auto face_id : cell.facesIds()) {
        const auto& face = mesh.face(face_id);

        // This is a boundary face
        if (face.isBoundary()) {
            grad += boundaryFaceIntegral(face);
            continue;
        }

        // This is an internal face
        // Area normal vector, poitning out of the cell
        auto Sf = mesh::outwardAreaVector(face, cell);
        const auto& nei = this->field()->mesh().otherSharingCell(cell, face);
        auto face_phi =
            0.5 * (this->field()->valueAtCell(cell) + this->field()->valueAtCell(nei));

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
    const auto& boundary_patch = this->field()->mesh().boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(this->field()->name());

    if (boundary_patch.isEmpty()) {
        return {0.0, 0.0, 0.0};
    }

    auto phi = this->field()->valueAtFace(face);
    return phi * face.area_vector();
}

LeastSquares::LeastSquares(field::IScalar* field) : IGradient(field) {
    const auto& mesh = this->field()->mesh();

    for (const auto& cell : mesh.cells()) {
        _cell_gradients.emplace_back(Vector3d::Zero());
    }
    setPseudoInvMatrices();
}
} // namespace prism::gradient