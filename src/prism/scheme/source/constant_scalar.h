#pragma once

#include "source.h"

namespace prism::scheme::source {
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

template <Sign SourceSign, field::IScalarBased Field>
ConstantScalar<SourceSign, Field>::ConstantScalar(Field phi)
    : _phi(phi), IExplicitSource(phi.mesh()->cellCount()) {}

template <Sign SourceSign, field::IScalarBased Field>
void ConstantScalar<SourceSign, Field>::apply() {
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = _phi.values().array() * vol_field.array();
        return;
    }
    rhs() = -_phi.values().array() * vol_field.array();
}
} // namespace prism::scheme::source
