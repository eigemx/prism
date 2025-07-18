#include "boundary.h"

namespace prism::gradient::boundary {
auto Symmetry::get(const prism::field::IScalar& field, const prism::mesh::Face& face) // NOLINT
    -> prism::Vector3d {
    return {0.0, 0.0, 0.0};
}

auto Outlet::get(const prism::field::IScalar& field, const prism::mesh::Face& face) // NOLINT
    -> prism::Vector3d {
    return {0.0, 0.0, 0.0};
}

auto Fixed::get(const prism::field::IScalar& field,
                const prism::mesh::Face& face) -> prism::Vector3d {
    /*
const auto& owner = field.mesh()->cell(face.owner());
prism::Vector3d d_Cf = face.center() - owner.center();
double d_Cf_norm = d_Cf.norm();
prism::Vector3d e = d_Cf / d_Cf_norm;

double delta_phi = field.valueAtFace(face) - field.valueAtCell(owner);
return (delta_phi / d_Cf_norm) * e;
*/
    const auto& owner = field.mesh()->cell(face.owner());
    const Vector3d grad_c = field.gradAtCellStored(owner);
    const Vector3d grad_b = grad_c - grad_c.dot(face.normal()) * face.normal();

    return grad_b;
}

auto NoSlip::get(const prism::field::IScalar& field,
                 const prism::mesh::Face& face) -> prism::Vector3d {
    Fixed fixed;
    return fixed.get(field, face);
}

auto VelocityInlet::get(const prism::field::IScalar& field,
                        const prism::mesh::Face& face) -> prism::Vector3d {
    Fixed fixed;
    return fixed.get(field, face);
}

auto FixedGradient::get(const prism::field::IScalar& field,
                        const prism::mesh::Face& face) -> prism::Vector3d {
    const auto& boundary_patch = field.mesh()->boundaryPatch(face);
    return boundary_patch.getVectorBoundaryCondition(field.name());
}

auto ZeroGradient::get(const prism::field::IScalar& field,                 // NOLINT
                       const prism::mesh::Face& face) -> prism::Vector3d { // NOLINT
    return {0.0, 0.0, 0.0};
}
} // namespace prism::gradient::boundary
