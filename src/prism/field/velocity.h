#pragma once

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
using Velocity = GeneralVector<VelocityComponent>;


} // namespace prism::field