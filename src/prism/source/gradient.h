#pragma once

#include <stdexcept>

#include "../field.h"
#include "../fvscheme.h"
#include "prism/gradient/gradient.h"
#include "prism/source/source.h"
#include "prism/types.h"

namespace prism::source {

template <SourceSign Sign, typename GradientScheme>
class Gradient : public FVScheme {
  public:
    Gradient(ScalarField& phi, Coords coord)
        : _phi(phi), _grad_scheme(phi), _coords(coord), FVScheme(phi.mesh().n_cells()) {}

    void apply() override;

    auto inline field() -> ScalarField& override { return _phi; }

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField& _phi;
    Coords _coords;
    GradientScheme _grad_scheme;
};

template <SourceSign Sign, typename GradientScheme>
void Gradient<Sign, GradientScheme>::apply() {
    auto grad_field = _grad_scheme.gradient_field();

    switch (_coords) {
        case Coords::X: {
            rhs() = grad_field.x().data();
            break;
        }
        case Coords::Y: {
            rhs() = grad_field.y().data();
            break;
        }

        case Coords::Z: {
            rhs() = grad_field.z().data();
            break;
        }

        default: {
            // We should not reach this!
            throw std::runtime_error(
                "source::Gradient::apply() was given an unknown coordinate!");
        }
    }

    if (Sign == SourceSign::Negative) {
        rhs() = -rhs();
    }
}


} // namespace prism::source