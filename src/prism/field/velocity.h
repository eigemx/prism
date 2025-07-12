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

/// TODO: for now VelocityBHManagerSetter is not needed, and we use general scalar field boundary
/// handlers manager setter. We can implement the below in case special treatment of velocity
/// fields at boundaries is needed.
class VelocityBHManagerSetter {
  public:
    using IVectorBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::vector::IVectorBoundaryHandler>;

    static void set(IVectorBHManager& manager);
};

using VelocityComponent = GeneralScalar<units::VelocityUnit, ScalarBHManagerSetter>;
using Velocity = GeneralVector<VelocityComponent, VectorBHManagerSetter>;


} // namespace prism::field