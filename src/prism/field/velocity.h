#pragma once

#include <array>
#include <string>

#include "scalar.h"
#include "units.h"
#include "vector.h"

namespace prism::field {

class VelocityCompBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using VelocityComponent = GeneralScalar<units::VelocityUnit, VelocityCompBHManagerSetter>;

class Velocity : public GeneralVector<VelocityComponent>, public units::VelocityUnit {
  public:
    Velocity(std::string name, const SharedPtr<mesh::PMesh>&, double value);
    Velocity(std::string name, const SharedPtr<mesh::PMesh>&, const Vector3d& data);
    Velocity(std::string name,
             const SharedPtr<mesh::PMesh>&,
             std::array<VelocityComponent, 3>& fields);

    using ComponentType = VelocityComponent;
};

} // namespace prism::field