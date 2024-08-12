#include "boundary.h"

#include <spdlog/spdlog.h>

#include "prism/mesh/face.h"
#include "scalar.h"

namespace prism::field::boundary {
auto Fixed::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh().boundary_patch(face);
    return patch.getScalarBoundaryCondition(field.name());
}

auto NoSlip::get(const field::Scalar& field, const mesh::Face& face) -> double {
    Fixed fixed;
    return fixed.get(field, face);
}

auto VelocityInlet::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh().boundary_patch(face);
    return patch.getScalarBoundaryCondition(field.name());
}

auto Empty::get(const field::Scalar& field, const mesh::Face& face) -> double {
    // TODO: Empty faces field value should not contribute to the solution, we need to avoid
    // having an "Empty" handler for fields.
    return field.values()[face.owner()];
}

auto Symmetry::get(const field::Scalar& field, const mesh::Face& face) -> double {
    return field.values()[face.owner()];
}

auto Outlet::get(const field::Scalar& field, const mesh::Face& face) -> double {
    return field.values()[face.owner()];
}

auto FixedGradient::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& mesh = field.mesh();
    const auto& patch = mesh.boundary_patch(face);

    Vector3d grad_at_boundary = patch.getVectorBoundaryCondition(field.name());
    const auto& owner = mesh.cell(face.owner());

    Vector3d e = face.center() - owner.center();
    double d_Cf = e.norm();
    e = e / e.norm();
    grad_at_boundary = grad_at_boundary * d_Cf;

    return grad_at_boundary.dot(e) + field.valueAtCell(owner);
}
} // namespace prism::field::boundary