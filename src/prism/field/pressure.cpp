#include "pressure.h"

#include "prism/log.h"

namespace prism::field {
void PressureBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::PressureBHManagerSetter(): setting default boundary handlers for a "
        "Pressure field");
    manager.addHandler<field::boundary::scalar::Fixed<Pressure>>();
    manager.addHandler<field::boundary::scalar::Symmetry<Pressure>>();
    manager.addHandler<field::boundary::scalar::Outlet<Pressure>>();
    manager.addHandler<field::boundary::scalar::ZeroGradient<Pressure>>();
}

} // namespace prism::field