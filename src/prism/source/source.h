#pragma once

#include <deque>

#include "../field.h"
#include "../fvscheme.h"

namespace prism::source {
// source term is assumed to be always on the right hand side of the conserved equation.

enum class SourceSign { Positive, Negative };

class ConstantScalar : public FVScheme {
  public:
    ConstantScalar(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}
    ConstantScalar(ScalarField&& phi) = delete;

    auto requires_correction() const -> bool override { return false; }

  private:
    // TODO: initializing rhs_vector() in constructor would suffice
    // NOLINTNEXTLINE (`face` is unused, but required by interface)
    void inline apply_interior(const mesh::Cell& cell, const mesh::Face& face) override {
        auto q = _phi[cell.id()];
        rhs_vector()(cell.id()) = q * cell.volume();
    }
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField& _phi;
};

template <SourceSign Sign = SourceSign::Positive>
class ImplicitPhi : public FVScheme {
  public:
    ImplicitPhi(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}
    ImplicitPhi(ScalarField&& phi) = delete;

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    auto inline requires_correction() const -> bool override { return false; }

    ScalarField& _phi;
};

} // namespace prism::source