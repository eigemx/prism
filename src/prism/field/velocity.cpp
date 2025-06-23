#include "velocity.h"

#include "prism/log.h"
#include "velocity_boundary.h"

namespace prism::field {
void VelocityCompBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::VelocityCompBHManagerSetter::set(): adding default boundary handlers for "
        "a VelocityComponent field instance");

    manager.addHandler<field::boundary::Fixed<VelocityComponent>>();
    manager.addHandler<field::boundary::VelocityInlet<VelocityComponent>>();
    manager.addHandler<field::boundary::Symmetry<VelocityComponent>>();
    manager.addHandler<field::boundary::Outlet<VelocityComponent>>();
    manager.addHandler<field::boundary::FixedGradient<VelocityComponent>>();
    manager.addHandler<field::boundary::NoSlip<VelocityComponent>>();
    manager.addHandler<field::boundary::ZeroGradient<VelocityComponent>>();
}

} // namespace prism::field
