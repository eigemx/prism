#pragma once

#include <optional>
#include <string>

#include "../field.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {
class GradientSchemeBase {
  public:
    virtual auto gradient(const mesh::Cell& c, const ScalarField& field) -> Vector3d = 0;
};

class GreenGauss : public GradientSchemeBase {
  public:
    auto gradient(const mesh::Cell& cell, const ScalarField& field) -> Vector3d override;
};
} // namespace prism::gradient