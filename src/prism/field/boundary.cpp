#include "boundary.h"

#include "field.h"
#include "prism/mesh/face.h"

namespace prism::field::boundary {
auto Fixed::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh().boundary_patch(face);
    return patch.get_scalar_bc(field.name());
}

auto VelocityInlet::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh().boundary_patch(face);
    return patch.get_scalar_bc(field.name());
}

auto Empty::get(const field::Scalar& field, const mesh::Face& face) -> double {
    return field.data()[face.owner()];
}

auto Symmetry::get(const field::Scalar& field, const mesh::Face& face) -> double {
    return field.data()[face.owner()];
}

auto Outlet::get(const field::Scalar& field, const mesh::Face& face) -> double {
    return field.data()[face.owner()];
}

auto FixedGradient::get(const field::Scalar& field, const mesh::Face& face) -> double {
    const auto& mesh = field.mesh();
    const auto& patch = mesh.boundary_patch(face);

    Vector3d grad_at_boundary = patch.get_vector_bc(name());
    const auto& owner = mesh.cell(face.owner());

    Vector3d e = face.center() - owner.center();
    double d_Cf = e.norm();
    e = e / e.norm();
    grad_at_boundary = grad_at_boundary * d_Cf;

    return grad_at_boundary.dot(e) + field.value_at_cell(owner);
}
} // namespace prism::field::boundary