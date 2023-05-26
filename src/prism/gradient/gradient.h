#pragma once

#include <vector>

#include "../field.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class GradientSchemeBase {
  public:
    GradientSchemeBase() = default;
    GradientSchemeBase(const GradientSchemeBase& g) = default;
    GradientSchemeBase(GradientSchemeBase&& g) = default;
    auto operator=(GradientSchemeBase&& g) -> GradientSchemeBase& = default;
    auto operator=(const GradientSchemeBase& g) -> GradientSchemeBase& = default;
    virtual ~GradientSchemeBase() = default;

    virtual auto gradient(const mesh::Cell& c) -> Vector3d = 0;
};

class GreenGauss : public GradientSchemeBase {
  public:
    GreenGauss(const ScalarField& field);
    auto gradient(const mesh::Cell& cell) -> Vector3d override;

  private:
    const ScalarField& _field;
    std::vector<Vector3d> _cell_gradients;
};


} // namespace prism::gradient