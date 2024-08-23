#include <cassert>

#include "pressure.h"

namespace prism::field::boundary {
template <>
auto NoSlip<Pressure>::get(const IScalar& field, const mesh::Face& face) -> double {
    assert(field.gradScheme() != nullptr);
    const auto& owner = field.mesh().cell(face.owner());
    const Vector3d d_Cb = face.center() - owner.center();

    // we use gradAtCellStored instead of gradAtCell, because gradAtCell invokes gradAtFace for
    // the boundary face at some point, which in turn will invoke NoSlip<Pressure>::get(), and
    // this will cause an infinite loop, to avoid this we use gradAtCellStored which uses the
    // gradient of the field at the cell calculated from the previous run of gradAtCell function.
    return field.valueAtCell(owner) + field.gradScheme()->gradAtCellStored(owner).dot(d_Cb);
}
} // namespace prism::field::boundary