#include "velocity_boundary.h"

namespace prism::field::boundary {

auto NoSlip<VelocityComponent>::get(const IScalar& field, const mesh::Face& face) -> double {
    Fixed<VelocityComponent> fixed;
    return fixed.get(field, face);
}

auto VelocityInlet<VelocityComponent>::get(const IScalar& field, const mesh::Face& face)
    -> double {
    Fixed<VelocityComponent> fixed;
    return fixed.get(field, face);
}
} // namespace prism::field::boundary