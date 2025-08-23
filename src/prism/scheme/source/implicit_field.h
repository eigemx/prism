#pragma once

#include "source.h"

namespace prism::scheme::source {
template <Sign SourceSign>
class ImplicitField : public IFullScheme, public IImplicitSource {
  public:
    ImplicitField(SharedPtr<field::Scalar>& phi) : IFullScheme(phi) {}
    ImplicitField(f64 coeff, SharedPtr<field::Scalar>& phi) : IFullScheme(phi), _coeff(coeff) {}

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    void applyInterior(const mesh::Face& face) override {}
    void applyBoundary() override {}

    f64 _coeff {1.0};
};


template <Sign SourceSign>
void ImplicitField<SourceSign>::apply() {
    matrix().setIdentity();
    matrix().diagonal() = _coeff * field()->mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        matrix() *= -1;
        return;
    }
}
} // namespace prism::scheme::source
