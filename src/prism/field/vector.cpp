#include "vector.h"

#include "prism/field/vector_boundary.h"

namespace prism::field {
void VectorBHManagerSetter::set(IVectorBHManager& manager) {
    log::debug(
        "prism::field::VectorBHManagerSetter::set(): adding default boundary handlers for a "
        "vector field instance");
    manager.addHandler<field::boundary::vector::Fixed<Vector>>();
    manager.addHandler<field::boundary::vector::NoSlip<Vector>>();
    manager.addHandler<field::boundary::vector::Symmetry<Vector>>();
    manager.addHandler<field::boundary::vector::Outlet<Vector>>();
    manager.addHandler<field::boundary::vector::ZeroGradient<Vector>>();
}
} // namespace prism::field