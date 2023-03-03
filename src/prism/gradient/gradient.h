#pragma once

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {
class GradientSchemeBase {};

class GreenGauss : public FVScheme, public GradientSchemeBase {
  public:
    GreenGauss() = default;

  private:
    auto apply_interior(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh) const -> AlteredCoeffs override;

    auto apply_boundary(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh) const -> AlteredCoeffs override;
};
} // namespace prism::gradient