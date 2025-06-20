#include "gradient.h"

#include "prism/mesh/utilities.h"

namespace prism::gradient {
IGradient::IGradient(field::IScalar* field) {
    if (field == nullptr) {
        throw std::runtime_error("prism::gradient::IGradient() was given a null IScalar pointer");
    }

    _field = field;
    log::debug("prism::gradient::IGradient() adding default boundary handlers for IGradient");
    this->boundaryHandlersManager().template addHandler<boundary::Fixed>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient>();
    this->boundaryHandlersManager().template addHandler<boundary::VelocityInlet>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip>();
}

auto IGradient::gradAtFace(const mesh::Face& face) -> Vector3d {
    // interpolate gradient at surrounding cells to the face center
    if (face.isInterior()) {
        // interior face
        const auto& mesh = _field->mesh();
        const auto& owner_cell = mesh->cell(face.owner());
        auto owner_grad = gradAtCell(owner_cell);

        const auto& neighbor_cell = mesh->cell(face.neighbor().value());
        auto neighbor_grad = gradAtCell(neighbor_cell);

        auto gc = mesh::geometricWeight(owner_cell, neighbor_cell, face);

        // Equation 9.33 without the correction part, a simple linear interpolation.
        return (gc * owner_grad) + ((1. - gc) * neighbor_grad);
    }

    // boundary face
    return gradAtBoundaryFace(face);
}

auto IGradient::gradAtBoundaryFace(const mesh::Face& face) -> Vector3d {
    const auto& boundary_patch = _field->mesh()->boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(_field->name());

    auto handler = this->boundaryHandlersManager().getHandler(boundary_condition.kindString());

    if (handler == nullptr) {
        throw prism::error::NonImplementedBoundaryCondition(
            "prism::gradient::IGradient::gradAtBoundaryFace()",
            boundary_patch.name(),
            boundary_condition.kindString());
    }

    return handler->get(*_field, face);
}
} // namespace prism::gradient