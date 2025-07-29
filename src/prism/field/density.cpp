#include "density.h"
#include "scalar_boundary.h"
#include "prism/log.h"

namespace prism::field {
void DensityBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::DensityBHManagerSetter(): setting default boundary handlers for a "
        "Density field");
    manager.addHandler<field::boundary::scalar::Fixed<Density>>();
    manager.addHandler<field::boundary::scalar::Symmetry<Density>>();
    manager.addHandler<field::boundary::scalar::Outlet<Density>>();
    manager.addHandler<field::boundary::scalar::ZeroGradient<Density>>();
}

} // namespace prism::field
