#include "green_gauss.h"

#include "prism/mesh/utilities.h"


namespace prism::gradient {
namespace {
auto boundaryFaceIntegral(const mesh::Face& face, const field::IScalar& field) -> Vector3d;
} // namespace

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


GreenGauss::GreenGauss(const SharedPtr<mesh::PMesh>& mesh) {
    _cell_gradients.resize(mesh->cellCount(), prism::Vector3d::Zero());
}

auto GreenGauss::gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d {
    // Lazy cache: recompute over the whole mesh only when the field's cell values changed.
    if (field.eventNo() != _computed_event) {
        computeAllCellGradients(field);
    }
    return _cell_gradients[cell.id()];
}

void GreenGauss::computeAllCellGradients(const field::IScalar& field) {
    for (const auto& cell : field.mesh()->cells()) {
        computeGradAtCell(cell, field);
    }
    _computed_event = field.eventNo();
}

auto GreenGauss::computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
    -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = field.mesh();
    const auto& phi = field;

    for (auto face_id : cell.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);

        // Empty faces contribute nothing to the gradient.
        if (face.isEmpty()) {
            continue;
        }

        if (face.isBoundary()) {
            grad += boundaryFaceIntegral(face, field);
            continue;
        }

        auto Sf = mesh::outwardAreaVector(face, cell);
        const auto& nei = field.mesh()->otherSharingCell(cell, face);
        auto face_phi = 0.5 * (phi.valueAtCell(cell) + phi.valueAtCell(nei));

        // skewness correction
        face_phi += correctSkewness(face, cell, nei);

        grad += Sf * face_phi;
    }
    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients[cell.id()] = grad;

    return grad;
}

namespace {
auto boundaryFaceIntegral(const mesh::Face& face, const field::IScalar& field) -> Vector3d {
    const auto& boundary_patch = field.mesh()->boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(field.name());
    auto phi = field.valueAtFace(face);
    return phi * face.areaVector();
}
} // namespace

} // namespace prism::gradient
