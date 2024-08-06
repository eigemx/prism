#pragma once

#include "fvscheme.h"
#include "prism/field/field.h"
#include "prism/gradient/gradient.h"
#include "prism/operations/operations.h"


namespace prism::scheme::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the SourceSign template parameter.
enum class SourceSign { Positive, Negative };

class ISource {};
class IImplicitSource : public ISource {};
class IExplicitSource : public ISource, public FVScheme<field::Scalar> {
  public:
    IExplicitSource(std::size_t n_cells);
    auto needsCorrection() const -> bool final { return false; }

  private:
    void inline apply_interior(const mesh::Face& face) final {}
    void inline apply_boundary(const mesh::Face& face) final {}
};

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
// coefficients are calculated during (and only during) initialization,
// no corrections are required afterwards

// TODO: ConstantScalar constructor should accept just a scalar value, and we should do the
// remaining housekeeping with creating the needed ScalarField
template <SourceSign Sign = SourceSign::Positive>
class ConstantScalar : public IExplicitSource {
  public:
    ConstantScalar(field::Scalar phi);
    void apply() override;

  private:
    field::Scalar _phi;
};

template <SourceSign Sign = SourceSign::Positive>
class Divergence : public IExplicitSource {
  public:
    Divergence(field::Vector& U);
    void apply() override;

  private:
    field::Vector _U;
};

// Adds a source for a gradient of a scalar field, in a specific coordinate
// for example the gradient of the pressure in the x-direction: ∂p/∂x
template <SourceSign Sign = SourceSign::Positive,
          typename GradientScheme = gradient::LeastSquares>
class Gradient : public IExplicitSource {
  public:
    Gradient(field::Scalar& phi, Coord coord);
    void apply() override;

  private:
    field::Scalar _phi;
    Coord _coord;
    GradientScheme _grad_scheme;
};

template <SourceSign Sign = SourceSign::Positive,
          typename GradientScheme = gradient::LeastSquares>
class Laplacian : public IExplicitSource {
  public:
    Laplacian(double kappa, field::Scalar phi);
    void apply() override;

  private:
    double _kappa;
    field::Scalar _phi;
    GradientScheme _grad_scheme;
};

// TODO: Test this!
template <SourceSign Sign = SourceSign::Positive>
class ImplicitField : public FVScheme<field::Scalar>, public IImplicitSource {
  public:
    ImplicitField(field::Scalar& phi) : _phi(phi), FVScheme(phi.mesh().nCells()) {}
    void apply() override;
    auto inline field() -> std::optional<field::Scalar> override { return _phi; }

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    auto inline needsCorrection() const -> bool override { return false; }

    field::Scalar _phi;
};

inline IExplicitSource::IExplicitSource(std::size_t n_cells) : FVScheme(n_cells, false) {}

template <SourceSign Sign>
ConstantScalar<Sign>::ConstantScalar(field::Scalar phi)
    : _phi(phi), IExplicitSource(phi.mesh().nCells()) {}

template <SourceSign Sign>
void inline ConstantScalar<Sign>::apply() {
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    if (Sign == SourceSign::Positive) {
        rhs() = _phi.values().array() * vol_field.array();
        return;
    }
    rhs() = -_phi.values().array() * vol_field.array();
}

template <SourceSign Sign>
Divergence<Sign>::Divergence(field::Vector& U) : IExplicitSource(U.mesh().nCells()), _U(U) {}

template <SourceSign Sign, typename GradientScheme>
Gradient<Sign, GradientScheme>::Gradient(field::Scalar& phi, Coord coord)
    : _phi(phi), _grad_scheme(phi), _coord(coord), IExplicitSource(phi.mesh().nCells()) {}

template <SourceSign Sign>
void inline Divergence<Sign>::apply() {
    if (Sign == SourceSign::Positive) {
        rhs() = ops::div(_U).values();
        return;
    }
    rhs() = -ops::div(_U).values();
}

template <SourceSign Sign, typename GradientScheme>
void Gradient<Sign, GradientScheme>::apply() {
    auto grad_field = _grad_scheme.gradient_field();
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    switch (_coord) {
        case Coord::X: {
            rhs() = grad_field.x().data().array() * vol_field.array();
            break;
        }
        case Coord::Y: {
            rhs() = grad_field.y().data().array() * vol_field.array();
            break;
        }

        case Coord::Z: {
            rhs() = grad_field.z().data().array() * vol_field.array();
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
Laplacian<Sign, GradientScheme>::Laplacian(double kappa, field::Scalar phi)
    : IExplicitSource(phi.mesh().nCells()),
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
void inline ImplicitField<Sign>::apply() {
    matrix().setIdentity();

    if (Sign == SourceSign::Positive) {
        matrix() *= -1;
        return;
    }
}

} // namespace prism::scheme::source