#pragma once

#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

class DiffusionSchemeBase {};

class Linear : public FVScheme, public DiffusionSchemeBase {
  public:
    Linear(double kappa, ScalarField& phi);
    inline void finalize() override {}
    auto mesh() const -> const mesh::PMesh& override { return _mesh; }


  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);

    double _kappa;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
};

} // namespace prism::diffusion