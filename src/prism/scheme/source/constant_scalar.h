#pragma once

#include "prism/field/scalar.h"
#include "source.h"

namespace prism::scheme::source {
// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
// coefficients are calculated during (and only during) initialization,
// no corrections are required afterwards
template <Sign SourceSign = Sign::Positive>
class ConstantScalar : public IExplicitSource {
  public:
    ConstantScalar(const SharedPtr<field::Scalar>& phi);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }
    auto field() const noexcept -> const SharedPtr<field::Scalar>& { return _phi; }

  private:
    SharedPtr<field::Scalar> _phi;
};

template <Sign SourceSign>
ConstantScalar<SourceSign>::ConstantScalar(const SharedPtr<field::Scalar>& phi)
    : _phi(phi), IExplicitSource(phi->mesh()->cellCount()) {}

template <Sign SourceSign>
void ConstantScalar<SourceSign>::apply() {
    const auto& vol_field = this->field()->mesh()->cellsVolumeVector();
    auto phi = this->field()->values();
    rhs() = phi.cwiseProduct(vol_field);

    if constexpr (SourceSign == Sign::Negative) {
        rhs() = -rhs();
        return;
    }
}

} // namespace prism::scheme::source
