#pragma once

#include "scalar.h"
#include "vector.h"

namespace prism::field {
/// TODO: for now VelocityBHManagerSetter is not needed, and we use general scalar field boundary
/// handlers manager setter. We can implement the below in case special treatment of velocity
/// fields at boundaries is needed.
class VelocityBHManagerSetter {
  public:
    using IVectorBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::vector::IVectorBoundaryHandler>;

    static void set(IVectorBHManager& manager);
};

class VelocityComponent : public Scalar {
  public:
    using Scalar::Scalar;
};

using Velocity = GeneralVector<VelocityComponent, VectorBHManagerSetter>;


} // namespace prism::field
