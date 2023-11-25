#pragma once

#include "fvscheme.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/operations/operations.h"


namespace prism::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the SourceSign template parameter.
enum class SourceSign { Positive, Negative };

class AbstractSource {};
class AbstractImplicitSource : public AbstractSource {};
class AbstractExplicitSource : public AbstractSource {};

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
// coefficients are calculated during (and only during) initialization,
// no corrections are required afterwards

// TODO: ConstantScalar constructor should accept just a scalar value, and we should do the
// remaining housekeeping with creating the needed ScalarField
template <SourceSign Sign = SourceSign::Positive>
class ConstantScalar : public FVScheme, public AbstractExplicitSource {
  public:
    ConstantScalar(ScalarField phi);

    auto requires_correction() const -> bool override { return false; }
    void apply() override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    const ScalarField _phi;
    VectorXd _volume_field;
};

template <SourceSign Sign = SourceSign::Positive>
class Divergence : public FVScheme, public AbstractExplicitSource {
  public:
    Divergence(const VectorField& U);

    void apply() override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    const VectorField _U;
};

// Adds a source for a gradient of a scalar field, in a specific coordinate
// for example the gradient of the pressure in the x-direction: ∂p/∂x
template <SourceSign Sign = SourceSign::Positive,
          typename GradientScheme = gradient::LeastSquares>
class Gradient : public FVScheme, public AbstractExplicitSource {
  public:
    Gradient(ScalarField& phi, Coord coord);

    void apply() override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    const ScalarField _phi;
    Coord _coord;
    GradientScheme _grad_scheme;
};

template <SourceSign Sign = SourceSign::Positive,
          typename GradientScheme = gradient::LeastSquares>
class Laplacian : public FVScheme, public AbstractExplicitSource {
  public:
    Laplacian(double kappa, ScalarField& phi);

    void apply() override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    double _kappa;
    const ScalarField _phi;
    GradientScheme _grad_scheme;
};

// TODO: Test this!
template <SourceSign Sign = SourceSign::Positive>
class Field : public FVScheme, public AbstractImplicitSource {
  public:
    Field(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {}

    void apply() override;
    auto inline field() -> std::optional<ScalarField> override { return _phi; }

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    auto inline requires_correction() const -> bool override { return false; }

    const ScalarField _phi;
};

template <SourceSign Sign>
ConstantScalar<Sign>::ConstantScalar(ScalarField phi)
    : _phi(phi), FVScheme(phi.mesh().n_cells(), false) {
    _volume_field.resize(phi.mesh().n_cells());

    for (const auto& cell : phi.mesh().cells()) {
        _volume_field[cell.id()] = cell.volume();
    }
}

template <SourceSign Sign>
void inline ConstantScalar<Sign>::apply() {
    if (Sign == SourceSign::Positive) {
        rhs() = _phi.data().array() * _volume_field.array();
        return;
    }

    rhs() = -_phi.data().array() * _volume_field.array();
}

template <SourceSign Sign>
Divergence<Sign>::Divergence(const VectorField& U) : FVScheme(U.mesh().n_cells(), false), _U(U) {}

template <SourceSign Sign, typename GradientScheme>
Gradient<Sign, GradientScheme>::Gradient(ScalarField& phi, Coord coord)
    : _phi(phi), _grad_scheme(phi), _coord(coord), FVScheme(phi.mesh().n_cells()) {}

template <SourceSign Sign>
void inline Divergence<Sign>::apply() {
    if (Sign == SourceSign::Positive) {
        rhs() = ops::div(_U).data();
        return;
    }
    rhs() = -ops::div(_U).data();
}

template <SourceSign Sign, typename GradientScheme>
void Gradient<Sign, GradientScheme>::apply() {
    auto grad_field = _grad_scheme.gradient_field();

    switch (_coord) {
        case Coord::X: {
            rhs() = grad_field.x().data();
            break;
        }
        case Coord::Y: {
            rhs() = grad_field.y().data();
            break;
        }

        case Coord::Z: {
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

template <SourceSign Sign, typename GradientScheme>
Laplacian<Sign, GradientScheme>::Laplacian(double kappa, ScalarField& phi)
    : FVScheme(phi.mesh().n_cells(), false),
      _kappa(kappa),
      _phi(phi),
      _grad_scheme(GradientScheme(phi)) {}

template <SourceSign Sign, typename GradientScheme>
void inline Laplacian<Sign, GradientScheme>::apply() {
    auto grad_phi = _grad_scheme.gradient_field();
    auto div = ops::div(grad_phi);

    if (Sign == SourceSign::Positive) {
        rhs() = div.data();
        return;
    }
    rhs() = -div.data();
}

template <SourceSign Sign>
void inline Field<Sign>::apply() {
    matrix().setIdentity();

    if (Sign == SourceSign::Positive) {
        matrix() *= -1;
        return;
    }
}

} // namespace prism::source