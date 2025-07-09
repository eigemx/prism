#include "velocity.h"

#include "prism/log.h"
#include "velocity_boundary.h"

namespace prism::field {
void VelocityCompBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::VelocityCompBHManagerSetter::set(): adding default boundary handlers for "
        "a VelocityComponent field instance");

    manager.addHandler<field::boundary::scalar::Fixed<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::VelocityInlet<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::Symmetry<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::Outlet<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::FixedGradient<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::NoSlip<VelocityComponent>>();
    manager.addHandler<field::boundary::scalar::ZeroGradient<VelocityComponent>>();
}

} // namespace prism::field
