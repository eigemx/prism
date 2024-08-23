#include "velocity.h"

#include "prism/log.h"
#include "velocity_boundary.h"

namespace prism::field {
Velocity::Velocity(std::string name, const mesh::PMesh& mesh, double value)
    : GeneralVector(std::move(name), mesh, value) {}

Velocity::Velocity(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : GeneralVector(std::move(name), mesh, data) {}

Velocity::Velocity(std::string name,
                   const mesh::PMesh& mesh,
                   std::array<VelocityComponent, 3>& fields)
    : GeneralVector(std::move(name), mesh, fields) {}

void VelocityCompBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::VelocityCompBHManagerSetter::set(): adding default boundary handlers for "
        "a VelocityComponent field instance");

    manager.addHandler<field::boundary::Fixed<VelocityComponent>>();
    manager.addHandler<field::boundary::VelocityInlet<VelocityComponent>>();
    manager.addHandler<field::boundary::Empty<VelocityComponent>>();
    manager.addHandler<field::boundary::Symmetry<VelocityComponent>>();
    manager.addHandler<field::boundary::Outlet<VelocityComponent>>();
    manager.addHandler<field::boundary::FixedGradient<VelocityComponent>>();
    manager.addHandler<field::boundary::NoSlip<VelocityComponent>>();
}

} // namespace prism::field