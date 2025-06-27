#include "pressure.h"

#include "prism/log.h"

namespace prism::field {
void PressureBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::PressureBHManagerSetter(): setting default boundary handlers for a "
        "Pressure field");
    manager.addHandler<field::boundary::Fixed<Pressure>>();
    manager.addHandler<field::boundary::Symmetry<Pressure>>();
    manager.addHandler<field::boundary::Outlet<Pressure>>();
    manager.addHandler<field::boundary::FixedGradient<Pressure>>();
    manager.addHandler<field::boundary::ZeroGradient<Pressure>>();
    // manager.addHandler<field::boundary::NoSlip<Pressure>>();
}

} // namespace prism::field