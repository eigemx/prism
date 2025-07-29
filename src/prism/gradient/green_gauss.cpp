#include "green_gauss.h"
#include "prism/mesh/utilities.h"

namespace prism::gradient {
auto GreenGauss::correctSkewness(const mesh::Face& face,
                                 const mesh::Cell& cell,
                                 const mesh::Cell& nei) const -> double {
    // This correction is based on option 1 from "Gradient Compuation" in Moukalled et al. (2016)
    // auto gc = mesh::geometricWeight(cell, nei, face);
    // double correction = gc * _cell_gradients[cell.id()].dot(face.center() - cell.center());
    // correction += (1.0 - gc) * _cell_gradients[nei.id()].dot(face.center() - nei.center());
    // return correction;

    // This correction is based on option 2 from "Gradient Compuation" in Moukalled et al. (2016)
    auto grad_sum = _cell_gradients[cell.id()] + _cell_gradients[nei.id()];
    auto vec = face.center() - (0.5 * (cell.center() + nei.center()));

    return 0.5 * grad_sum.dot(vec);
}


GreenGauss::GreenGauss(field::IScalar* field) : IGradient(field) {
    const auto& mesh = field->mesh();
    _cell_gradients.resize(mesh->cellCount(), prism::Vector3d::Zero());
}

auto GreenGauss ::gradAtCellStored(const mesh::Cell& cell) -> Vector3d {
    return _cell_gradients[cell.id()];
}

auto GreenGauss ::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = this->field()->mesh();
    const auto& phi = this->field();

    for (auto face_id : cell.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);

        if (face.isBoundary()) {
            grad += boundaryFaceIntegral(face);
            continue;
        }

        auto Sf = mesh::outwardAreaVector(face, cell);
        const auto& nei = this->field()->mesh()->otherSharingCell(cell, face);
        auto face_phi = 0.5 * (phi->valueAtCell(cell) + phi->valueAtCell(nei));

        // skewness correction
        face_phi += correctSkewness(face, cell, nei);

        grad += Sf * face_phi;
    }
    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients[cell.id()] = grad;

    return grad;
}

auto GreenGauss::boundaryFaceIntegral(const mesh::Face& face) -> Vector3d {
    const auto& boundary_patch = this->field()->mesh()->boundaryPatch(face);
    if (boundary_patch.isEmpty()) {
        return {0.0, 0.0, 0.0};
    }

    const auto& boundary_condition = boundary_patch.getBoundaryCondition(this->field()->name());
    auto phi = this->field()->valueAtFace(face);
    return phi * face.areaVector();
}

} // namespace prism::gradient
