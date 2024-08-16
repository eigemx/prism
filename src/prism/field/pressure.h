#pragma once

#include "prism/field/scalar.h"

namespace prism::field {

class PressureBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Pressure = GeneralScalar<units::Pascal, PressureBHManagerSetter>;

} // namespace prism::field