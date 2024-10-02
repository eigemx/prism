#pragma once

#include "prism/field/scalar.h"
#include "prism/log.h"
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
template <SourceSign Sign = SourceSign::Positive, typename Field = field::Scalar>
class Gradient : public IExplicitSource {
  public:
    Gradient(Field phi, Coord coord);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Field _phi;
    Coord _coord;
};

template <SourceSign Sign = SourceSign::Positive,
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

// TODO: Test this!
template <SourceSign Sign, field::IScalarBased Field>
class ImplicitField : public IFullScheme<Field>, public IImplicitSource {
  public:
    ImplicitField(Field& phi) : _phi(phi), IFullScheme<Field>(phi.mesh().nCells()) {}
    void apply() override;
    auto inline field() -> Field override { return _phi; }
    auto inline needsCorrection() const -> bool override { return false; }

  private:
    void inline applyInterior(const mesh::Face& face) override {}

    field::Scalar _phi;
};

template <SourceSign Sign, field::IScalarBased Field>
ConstantScalar<Sign, Field>::ConstantScalar(Field phi)
    : _phi(phi), IExplicitSource(phi.mesh().cellCount()) {}

template <SourceSign Sign, field::IScalarBased Field>
void inline ConstantScalar<Sign, Field>::apply() {
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    if constexpr (Sign == SourceSign::Positive) {
        rhs() = _phi.values().array() * vol_field.array();
        return;
    }
    rhs() = -_phi.values().array() * vol_field.array();
}

template <SourceSign Sign, typename Vector>
Divergence<Sign, Vector>::Divergence(Vector U) : IExplicitSource(U.mesh().cellCount()), _U(U) {}


template <SourceSign Sign, typename Vector>
void inline Divergence<Sign, Vector>::apply() {
    const auto& vol_field = _U.mesh().cellsVolumeVector();
    if constexpr (Sign == SourceSign::Positive) {
        rhs() = ops::div(_U).values().array() * vol_field.array();
        return;
    }
    rhs() = -ops::div(_U).values().array() * vol_field.array();
}

template <SourceSign Sign, typename Field>
Gradient<Sign, Field>::Gradient(Field phi, Coord coord)
    : _phi(phi), _coord(coord), IExplicitSource(phi.mesh().cellCount()) {
    log::debug("prism::scheme::source::Gradient(): Creating gradient source for field '{}'",
               phi.name());
}

template <SourceSign Sign, typename Field>
void Gradient<Sign, Field>::apply() {
    const auto& vol_field = _phi.mesh().cellsVolumeVector();

    switch (_coord) {
        case Coord::X: {
            rhs() = ops::grad(_phi, Coord::X).values().array() * vol_field.array();
            break;
        }
        case Coord::Y: {
            rhs() = ops::grad(_phi, Coord::Y).values().array() * vol_field.array();
            break;
        }

        case Coord::Z: {
            rhs() = ops::grad(_phi, Coord::Z).values().array() * vol_field.array();
            break;
        }
    }

    if constexpr (Sign == SourceSign::Negative) {
        rhs() = -rhs();
    }
}

template <SourceSign Sign, typename Kappa, typename Field>
Laplacian<Sign, Kappa, Field>::Laplacian(Kappa kappa, Field phi)
    : IExplicitSource(phi.mesh().nCells()), _kappa(kappa), _phi(phi) {}

template <SourceSign Sign, typename Kappa, typename Field>
void inline Laplacian<Sign, Kappa, Field>::apply() {
    auto grad_phi = ops::grad(_phi);
    auto div = ops::div(grad_phi);

    if constexpr (Sign == SourceSign::Positive) {
        rhs() = div.values();
        return;
    }
    rhs() = -div.values();
}

template <SourceSign Sign, field::IScalarBased Field>
void inline ImplicitField<Sign, Field>::apply() {
    this->matrix().setIdentity();

    if constexpr (Sign == SourceSign::Positive) {
        this->matrix() *= -1;
        return;
    }
}

} // namespace prism::scheme::source
