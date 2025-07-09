#pragma once

#include "scalar.h"
#include "units.h"
#include "vector.h"

namespace prism::field {

class VelocityCompBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::scalar::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

class VelocityBHManagerSetter {
  public:
    using IVectorBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::vector::IVectorBoundaryHandler>;

    static void set(IVectorBHManager& manager);
};

using VelocityComponent = GeneralScalar<units::VelocityUnit, VelocityCompBHManagerSetter>;
using Velocity = GeneralVector<VelocityComponent, VectorBHManagerSetter>;


} // namespace prism::field