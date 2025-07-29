#pragma once

#include "prism/field/scalar.h"

namespace prism::field {

class DensityBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::scalar::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Density = GeneralScalar<units::DensityUnit, DensityBHManagerSetter>;

} // namespace prism::field
