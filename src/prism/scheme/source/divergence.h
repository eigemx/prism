#pragma once

#include "prism/operations/operations.h"
#include "source.h"

namespace prism::scheme::source {

template <Sign SourceSign = Sign::Positive>
class Divergence : public IExplicitSource {
  public:
    Divergence(const SharedPtr<field::IVector>& U);

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    SharedPtr<field::IVector> _U;
};

template <Sign SourceSign>
Divergence<SourceSign>::Divergence(const SharedPtr<field::IVector>& U)
    : IExplicitSource(U->mesh()->cellCount()), _U(U) {}


template <Sign SourceSign>
void Divergence<SourceSign>::apply() {
    const auto& vol_field = _U->mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = ops::div(*_U).values().cwiseProduct(vol_field);
        return;
    }
    rhs() = -ops::div(*_U).values().cwiseProduct(vol_field);
}


} // namespace prism::scheme::source
