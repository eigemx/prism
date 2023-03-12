#pragma once

#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {
class GradientSchemeBase {};

class GreenGauss : public GradientSchemeBase {
  public:
    GreenGauss() = default;
};
} // namespace prism::gradient