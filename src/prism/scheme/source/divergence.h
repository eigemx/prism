#pragma once

#include "prism/field/vector.h"
#include "prism/operations/operations.h"
#include "source.h"

namespace prism::scheme::source {

template <Sign SourceSign = Sign::Positive, typename Vector = field::Vector>
class Divergence : public IExplicitSource {
  public:
    Divergence(Vector U);

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Vector _U;
};

template <Sign SourceSign, typename Vector>
Divergence<SourceSign, Vector>::Divergence(Vector U)
    : IExplicitSource(U.mesh()->cellCount()), _U(U) {}


template <Sign SourceSign, typename Vector>
void Divergence<SourceSign, Vector>::apply() {
    const auto& vol_field = _U.mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = ops::div(_U).values().array() * vol_field.array();
        return;
    }
    rhs() = -ops::div(_U).values().array() * vol_field.array();
}


} // namespace prism::scheme::source
