#include <cassert>

#include "pressure.h"

namespace prism::field::boundary::scalar {
template <>
auto Symmetry<Pressure>::get(const IScalar& field, const mesh::Face& face) -> double {
    // This function is based on equations (15.151), (15.152) and (15.153) from Moukallad et. al.
    // gradient in normal direction is zero, ∇pb.n = 0, so the value of the field should be
    // extrapolated. To make sure that we have a zero normal gradient, the field gradient at the
    // boundary is computed as:
    // ∇pb = ∇pc - (∇pc.n)n
    const auto& owner = field.mesh()->cell(face.owner());
    const Vector3d grad_c = field.gradAtCellStored(owner);
    const Vector3d grad_b = grad_c - grad_c.dot(face.normal()) * face.normal();

    // pb = pc + ∇pb.dCb
    const Vector3d dCb = face.center() - owner.center();
    return field.valueAtCell(face.owner()) + grad_b.dot(dCb);
}
} // namespace prism::field::boundary::scalar