#pragma once

#include "prism/field/ifield.h"
#include "prism/field/scalar.h"
#include "prism/log.h"
#include "prism/operations/operations.h"
#include "scheme.h"


namespace prism::scheme::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the Sign template parameter.

class ISource {};

class IImplicitSource : public ISource {};

template <typename T>
concept IImplicitSourceBased = std::derived_from<T, IImplicitSource>;

// We inherit from FVScheme with field::Scalar as template specialization because for explicit
// sources the type of the field won't matter, because we are not contributing to the matrix of
// coefficients for the linear system of the conserved equation.
class IExplicitSource : public ISource, public IPartialScheme {
  public:
    IExplicitSource(std::size_t n_cells) : IPartialScheme(n_cells) {}
};

template <typename T>
concept IExplicitSourceBased = std::derived_from<T, IExplicitSource>;

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
// coefficients are calculated during (and only during) initialization,
// no corrections are required afterwards

/// TODO: ConstantScalar constructor should accept just a scalar value, and we should do the
// remaining housekeeping with creating the needed ScalarField
template <Sign SourceSign = Sign::Positive, field::IScalarBased Field = field::Scalar>
class ConstantScalar : public IExplicitSource {
  public:
    ConstantScalar(Field phi);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    Field _phi;
};

template <Sign SourceSign = Sign::Positive, typename Vector = field::Vector>
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
template <Sign SourceSign = Sign::Positive, typename Field = field::Scalar>
class Gradient : public IExplicitSource {
  public:
    Gradient(Field phi, Coord coord);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Field _phi;
    Coord _coord;
};

template <Sign SourceSign = Sign::Positive,
          typename Kappa = field::UniformScalar,
          typename Field = field::Scalar>
class Laplacian : public IExplicitSource {
  public:
    Laplacian(Kappa kappa, Field phi);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Kappa _kappa;
    Field _phi;
};

template <Sign SourceSign, field::IScalarBased Field>
class ImplicitField : public IFullScheme<Field>, public IImplicitSource {
  public:
    ImplicitField(Field& phi) : _phi(phi), IFullScheme<Field>(phi.mesh()->cellCount()) {}
    ImplicitField(double coeff, Field& phi)
        : _phi(phi), IFullScheme<Field>(phi.mesh()->cellCount()), _coeff(coeff) {}

    void apply() override;
    auto inline field() -> Field override { return _phi; }
    auto inline needsCorrection() const noexcept -> bool override { return false; }

  private:
    void inline applyInterior(const mesh::Face& face) override {}
    void inline applyBoundary() override {}

    field::Scalar _phi;
    double _coeff {1.0};
};

template <Sign SourceSign, field::IScalarBased Field>
ConstantScalar<SourceSign, Field>::ConstantScalar(Field phi)
    : _phi(phi), IExplicitSource(phi.mesh()->cellCount()) {}

template <Sign SourceSign, field::IScalarBased Field>
void inline ConstantScalar<SourceSign, Field>::apply() {
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = _phi.values().array() * vol_field.array();
        return;
    }
    rhs() = -_phi.values().array() * vol_field.array();
}

template <Sign SourceSign, typename Vector>
Divergence<SourceSign, Vector>::Divergence(Vector U)
    : IExplicitSource(U.mesh()->cellCount()), _U(U) {}


template <Sign SourceSign, typename Vector>
void inline Divergence<SourceSign, Vector>::apply() {
    const auto& vol_field = _U.mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = ops::div(_U).values().array() * vol_field.array();
        return;
    }
    rhs() = -ops::div(_U).values().array() * vol_field.array();
}

template <Sign SourceSign, typename Field>
Gradient<SourceSign, Field>::Gradient(Field phi, Coord coord)
    : _phi(phi), _coord(coord), IExplicitSource(phi.mesh()->cellCount()) {
    log::debug(
        "prism::scheme::source::Gradient(): Creating {}-coordinate gradient source for field "
        "'{}'",
        field::coordToStr(coord),
        phi.name());
}

template <Sign SourceSign, typename Field>
void Gradient<SourceSign, Field>::apply() {
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();
    rhs() = ops::grad(_phi, _coord).values().cwiseProduct(vol_field);

    if constexpr (SourceSign == Sign::Negative) {
        rhs() = -rhs();
    }
}

template <Sign SourceSign, typename Kappa, typename Field>
Laplacian<SourceSign, Kappa, Field>::Laplacian(Kappa kappa, Field phi)
    : IExplicitSource(phi.mesh()->nCells()), _kappa(kappa), _phi(phi) {}

template <Sign SourceSign, typename Kappa, typename Field>
void inline Laplacian<SourceSign, Kappa, Field>::apply() {
    auto grad_phi = ops::grad(_phi);
    auto div = ops::div(grad_phi);

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = div.values();
        return;
    }
    rhs() = -div.values();
}

template <Sign SourceSign, field::IScalarBased Field>
void inline ImplicitField<SourceSign, Field>::apply() {
    this->matrix().setIdentity();
    this->matrix().diagonal() *= _coeff * this->field().mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        this->matrix() *= -1;
        return;
    }
}

} // namespace prism::scheme::source
