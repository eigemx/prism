#pragma once

#include "prism/field/scalar.h"
#include "prism/gradient/gradient.h"
#include "prism/operations/operations.h"
#include "scheme.h"


namespace prism::scheme::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the SourceSign template parameter.
enum class SourceSign { Positive, Negative };

class ISource {};

class IImplicitSource : public ISource {};

// We inherit from FVScheme with field::Scalar as template specialization because for explicit
// sources the type of the field won't matter, because we are not contributing to the matrix of
// coefficients for the linear system of the conserved equation.
class IExplicitSource : public ISource, public IPartialScheme {
  public:
    IExplicitSource(std::size_t n_cells) : IPartialScheme(n_cells) {}
};

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
// coefficients are calculated during (and only during) initialization,
// no corrections are required afterwards

// TODO: ConstantScalar constructor should accept just a scalar value, and we should do the
// remaining housekeeping with creating the needed ScalarField
template <SourceSign Sign = SourceSign::Positive, field::IScalarBased Field = field::Scalar>
class ConstantScalar : public IExplicitSource {
  public:
    ConstantScalar(Field phi);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    Field _phi;
};

template <SourceSign Sign = SourceSign::Positive, typename Vector = field::Vector>
class Divergence : public IExplicitSource {
  public:
    Divergence(Vector U);

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Vector _U;
};

// Adds a source for a gradient of a scalar field, in a specific coordinate
// for example the gradient of the pressure in the x-direction: ∂p/∂x
template <SourceSign Sign = SourceSign::Positive,
          typename Field = field::Scalar,
          typename GradientScheme = gradient::LeastSquares<Field>>
class Gradient : public IExplicitSource {
  public:
    Gradient(Field phi, Coord coord);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }


  private:
    Field _phi;
    Coord _coord;
    GradientScheme _grad_scheme;
};

template <SourceSign Sign = SourceSign::Positive,
          typename Kappa = field::UniformScalar,
          typename Field = field::Scalar,
          typename GradientScheme = gradient::LeastSquares<Field>>
class Laplacian : public IExplicitSource {
  public:
    Laplacian(Kappa kappa, Field phi);
    void apply() override;

  private:
    Kappa _kappa;
    Field _phi;
    GradientScheme _grad_scheme;
};

// TODO: Test this!
template <SourceSign Sign, field::IScalarBased Field>
class ImplicitField : public IFullScheme<Field>, public IImplicitSource {
  public:
    ImplicitField(Field& phi) : _phi(phi), IFullScheme<Field>(phi.mesh().nCells()) {}
    void apply() override;
    auto inline field() -> Field override { return _phi; }

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    auto inline needsCorrection() const -> bool override { return false; }

    field::Scalar _phi;
};

template <SourceSign Sign, field::IScalarBased Field>
ConstantScalar<Sign, Field>::ConstantScalar(Field phi)
    : _phi(phi), IExplicitSource(phi.mesh().nCells()) {}

template <SourceSign Sign, field::IScalarBased Field>
void inline ConstantScalar<Sign, Field>::apply() {
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    if (Sign == SourceSign::Positive) {
        rhs() = _phi.values().array() * vol_field.array();
        return;
    }
    rhs() = -_phi.values().array() * vol_field.array();
}

template <SourceSign Sign, typename Vector>
Divergence<Sign, Vector>::Divergence(Vector U) : IExplicitSource(U.mesh().nCells()), _U(U) {}


template <SourceSign Sign, typename Vector>
void inline Divergence<Sign, Vector>::apply() {
    if (Sign == SourceSign::Positive) {
        rhs() = ops::div(_U).values();
        return;
    }
    rhs() = -ops::div(_U).values();
}

template <SourceSign Sign, typename Field, typename GradientScheme>
Gradient<Sign, Field, GradientScheme>::Gradient(Field phi, Coord coord)
    : _phi(phi), _grad_scheme(phi), _coord(coord), IExplicitSource(phi.mesh().nCells()) {}

template <SourceSign Sign, typename Field, typename GradientScheme>
void Gradient<Sign, Field, GradientScheme>::apply() {
    auto grad_field = _grad_scheme.gradField();
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    switch (_coord) {
        case Coord::X: {
            rhs() = grad_field.x().values().array() * vol_field.array();
            break;
        }
        case Coord::Y: {
            rhs() = grad_field.y().values().array() * vol_field.array();
            break;
        }

        case Coord::Z: {
            rhs() = grad_field.z().values().array() * vol_field.array();
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

template <SourceSign Sign, typename Kappa, typename Field, typename GradientScheme>
Laplacian<Sign, Kappa, Field, GradientScheme>::Laplacian(Kappa kappa, Field phi)
    : IExplicitSource(phi.mesh().nCells()),
      _kappa(kappa),
      _phi(phi),
      _grad_scheme(GradientScheme(phi)) {}

template <SourceSign Sign, typename Kappa, typename Field, typename GradientScheme>
void inline Laplacian<Sign, Kappa, Field, GradientScheme>::apply() {
    auto grad_phi = _grad_scheme.gradField();
    auto div = ops::div(grad_phi);

    if (Sign == SourceSign::Positive) {
        rhs() = div.values();
        return;
    }
    rhs() = -div.values();
}

template <SourceSign Sign, field::IScalarBased Field>
void inline ImplicitField<Sign, Field>::apply() {
    this->matrix().setIdentity();

    if (Sign == SourceSign::Positive) {
        this->matrix() *= -1;
        return;
    }
}

} // namespace prism::scheme::source