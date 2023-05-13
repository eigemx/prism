#pragma once

#include <optional>
#include <string>
#include <vector>

#include "../field.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {
class GradientSchemeBase {
  public:
    GradientSchemeBase() = default;
    GradientSchemeBase(const GradientSchemeBase& g) = delete;
    GradientSchemeBase(GradientSchemeBase&& g) = delete;
    auto operator=(GradientSchemeBase&& g) -> GradientSchemeBase& = delete;
    auto operator=(const GradientSchemeBase& g) -> GradientSchemeBase& = delete;

    virtual ~GradientSchemeBase() = default;

    virtual auto gradient(const mesh::Cell& c) -> Vector3d = 0;
};

class GreenGauss : public GradientSchemeBase {
  public:
    GreenGauss(const ScalarField& field);
    auto gradient(const mesh::Cell& cell) -> Vector3d override;
    auto field_gradients() const -> const std::vector<Vector3d>& { return _cell_gradients; }
    auto field_gradients() -> std::vector<Vector3d>& { return _cell_gradients; }

  private:
    const ScalarField& _field;
    std::vector<Vector3d> _cell_gradients;
};
} // namespace prism::gradient