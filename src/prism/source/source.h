#pragma once

#include <deque>

#include "../field.h"
#include "../fvscheme.h"

namespace prism::source {
// General note: source term is assumed to be always on the right hand side of the conserved equation.

class ConstantScalar : public FVScheme {
  public:
    ConstantScalar(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}
    ConstantScalar(ScalarField&& phi) : _phi(phi), FVScheme(_phi.mesh().n_cells()) {}

    auto requires_correction() const -> bool override { return false; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField& _phi;
};
} // namespace prism::source