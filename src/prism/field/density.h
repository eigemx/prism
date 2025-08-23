#pragma once

#include "prism/field/scalar.h"

namespace prism::field {
class Density : public Scalar {
  public:
    using Scalar::Scalar;   // Inherit constructors of field::Scalar
};

} // namespace prism::field
