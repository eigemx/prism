#pragma once

#include <deque>

#include "../field.h"
#include "../fvscheme.h"

namespace prism::source {
// source term is assumed to be always on the right hand side of the conserved equation.

enum class SourceSign { Positive, Negative };

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
class ConstantScalar : public FVScheme {
  public:
    ConstantScalar(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}
    ConstantScalar(ScalarField&& phi) = delete;

    auto requires_correction() const -> bool override { return false; }

    auto inline field() -> ScalarField& override { return _phi; }

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

// TODO: Test this!
template <SourceSign Sign = SourceSign::Positive>
class ImplicitPhi : public FVScheme {
  public:
    ImplicitPhi(ScalarField& phi, double coeff = 1.0)
        : _phi(phi), _coeff(coeff), FVScheme(phi.mesh().n_cells()) {}
    ImplicitPhi(ScalarField&& phi) = delete;

    auto inline field() -> ScalarField& override { return _phi; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    auto inline requires_correction() const -> bool override { return false; }

    ScalarField& _phi;
    double _coeff {1.0};
};

} // namespace prism::source