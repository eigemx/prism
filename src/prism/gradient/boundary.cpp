#include "boundary.h"

namespace prism::gradient::boundary {
auto Symmetry::get(field::IScalar& field, const prism::mesh::Face& face) // NOLINT
    -> Vector3d {
    return field.gradAtCell(field.mesh()->cell(face.owner()));
}

auto Outlet::get(field::IScalar& field, const prism::mesh::Face& face) // NOLINT
    -> Vector3d {
    return field.gradAtCell(field.mesh()->cell(face.owner()));
}

auto Fixed::get(field::IScalar& field, const prism::mesh::Face& face) -> Vector3d {
    // The following is based on uFVM gradient correction in cfdUpdateGradient function
    const auto& owner = field.mesh()->cell(face.owner());
    const Vector3d grad_c = field.gradAtCell(owner);
    const Vector3d d_Cf = face.center() - owner.center();
    const Vector3d e = d_Cf / d_Cf.norm();

    double phi_b = field.valueAtFace(face);
    double phi_C = field.valueAtCell(owner);
    return grad_c - (grad_c.dot(e)) * e + (phi_b - phi_C) / d_Cf.norm() * e;
}

auto NoSlip::get(field::IScalar& field, const prism::mesh::Face& face) -> Vector3d {
    Fixed fixed;
    return fixed.get(field, face);
}

auto ZeroGradient::get(field::IScalar& field,                       // NOLINT
                       const prism::mesh::Face& face) -> Vector3d { // NOLINT
    return field.gradAtCell(field.mesh()->cell(face.owner()));
}
} // namespace prism::gradient::boundary
