#include "boundary.h"

namespace prism::gradient::boundary {
auto Empty::get(const prism::field::Scalar& field, const prism::mesh::Face& face) // NOLINT
    -> prism::Vector3d {
    return {0.0, 0.0, 0.0};
}

auto Symmetry::get(const prism::field::Scalar& field, const prism::mesh::Face& face) // NOLINT
    -> prism::Vector3d {
    return {0.0, 0.0, 0.0};
}

auto Outlet::get(const prism::field::Scalar& field, const prism::mesh::Face& face) // NOLINT
    -> prism::Vector3d {
    return {0.0, 0.0, 0.0};
}

auto Fixed::get(const prism::field::Scalar& field, const prism::mesh::Face& face)
    -> prism::Vector3d {
    const auto& owner = field.mesh().cell(face.owner());
    prism::Vector3d d_Cf = face.center() - owner.center();
    double d_Cf_norm = d_Cf.norm();
    prism::Vector3d e = d_Cf / d_Cf_norm;

    double delta_phi = field.value_at_face(face) - field.value_at_cell(owner);
    return (delta_phi / d_Cf_norm) * e;
}

auto VelocityInlet::get(const prism::field::Scalar& field, const prism::mesh::Face& face)
    -> prism::Vector3d {
    Fixed fixed;
    return fixed.get(field, face);
}

auto FixedGradient::get(const prism::field::Scalar& field, const prism::mesh::Face& face)
    -> prism::Vector3d {
    const auto& boundary_patch = field.mesh().boundary_patch(face);
    return boundary_patch.get_vector_bc(field.name());
}
} // namespace prism::gradient::boundary