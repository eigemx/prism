#include "green_gauss.h"

#include <algorithm>


namespace prism::gradient {
namespace {
auto boundaryFaceIntegral(const mesh::Face& face, const field::IScalar& field) -> Vector3d;
} // namespace

auto GreenGauss::skewnessCorrection(const mesh::Face& /*face*/,
                                    const mesh::Cell& /*cell*/,
                                    const mesh::Cell& /*nei*/) const -> f64 {
    return 0.0;
}


GreenGauss::GreenGauss(const SharedPtr<mesh::PMesh>& mesh) {
    _cell_gradients.resize(mesh->cellCount(), prism::Vector3d::Zero());
    _previous_gradients.resize(mesh->cellCount(), prism::Vector3d::Zero());
}

auto GreenGauss::gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d {
    // Lazy cache: recompute over the whole mesh only when the field's cell values changed.
    if (field.eventNo() != _computed_event) {
        computeAllCellGradients(field);
    }
    return _cell_gradients[cell.id()];
}

void GreenGauss::computeAllCellGradients(const field::IScalar& field) {
    const auto& mesh = field.mesh();
    const auto& volumes = mesh->cellsVolumeVector();

    // Keep the previous gradient field for the skewness correction, and reuse the existing
    // storage as the new gradient buffer.
    _previous_gradients.swap(_cell_gradients);
    std::fill(_cell_gradients.begin(), _cell_gradients.end(), Vector3d::Zero());

    for (const auto& face : mesh->interiorFaces()) {
        const auto& owner = mesh->cell(face.owner());
        const auto& neighbor = mesh->cell(face.neighbor().value());
        auto face_phi = 0.5 * (field.valueAtCell(owner) + field.valueAtCell(neighbor));

        // skewness correction
        face_phi += skewnessCorrection(face, owner, neighbor);

        const Vector3d contribution = face.areaVector() * face_phi;
        _cell_gradients[owner.id()] += contribution;
        _cell_gradients[neighbor.id()] -= contribution;
    }

    for (const auto& face : mesh->nonEmptyBoundaryFaces()) {
        _cell_gradients[face.owner()] += boundaryFaceIntegral(face, field);
    }

    for (size_t i = 0; i < _cell_gradients.size(); ++i) {
        _cell_gradients[i] /= volumes[static_cast<Eigen::Index>(i)];
    }

    _computed_event = field.eventNo();
}

auto GreenGauss::computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
    -> Vector3d {
    return gradAtCell(cell, field);
}

namespace {
auto boundaryFaceIntegral(const mesh::Face& face, const field::IScalar& field) -> Vector3d {
    const auto& boundary_patch = field.mesh()->boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(field.name());
    auto phi = field.valueAtFace(face);
    return phi * face.areaVector();
}
} // namespace

auto SkewCorrectedGreenGauss::skewnessCorrection(const mesh::Face& face,
                                                 const mesh::Cell& cell,
                                                 const mesh::Cell& nei) const -> f64 {
    // This correction is based on option 2 from "Gradient Computation" in Moukalled et al.
    // (2016).
    auto grad_sum = _previous_gradients[cell.id()] + _previous_gradients[nei.id()];
    auto vec = face.center() - (0.5 * (cell.center() + nei.center()));

    return 0.5 * grad_sum.dot(vec);
}

} // namespace prism::gradient
