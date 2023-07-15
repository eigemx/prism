#pragma once

#include <deque>

#include "../field.h"
#include "../fvscheme.h"

namespace prism::source {
// General note: source term is assumed to be always on the right hand side of the conserved equation.

class ConstantScalar : public FVScheme {
  public:
    ConstantScalar() = delete;
    ConstantScalar(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}
    ConstantScalar(ScalarField&& phi) = delete;

    ConstantScalar(const ConstantScalar& other) = delete;
    ConstantScalar(ConstantScalar&& other) = delete;
    auto operator=(const ConstantScalar& other) -> ConstantScalar& = delete;
    auto operator=(ConstantScalar&& other) -> ConstantScalar& = delete;
    ~ConstantScalar() = default;

    auto requires_correction() const -> bool override { return false; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField& _phi;
};

class Phi : public FVScheme {
  public:
    Phi(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()), _phi_old(phi) {}
    Phi(ScalarField&& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()), _phi_old(phi) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField& _phi;
    ScalarField _phi_old;
};

} // namespace prism::source