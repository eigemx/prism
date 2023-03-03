#pragma once

#include <memory>

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

class DiffusionSchemeBase {};

class Linear : public FVScheme, public DiffusionSchemeBase {
  public:
    Linear(double kappa) : _kappa(kappa) {}


  private:
    auto apply_interior(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh) const -> AlteredCoeffs override;

    auto apply_boundary(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh) const -> AlteredCoeffs override;

    auto apply_wall_boundary(const mesh::Cell& cell,
                             const mesh::Face& face,
                             double phi_wall) const -> AlteredCoeffs;

    double _kappa;
};

} // namespace prism::diffusion