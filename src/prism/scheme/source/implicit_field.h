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
    resetValues();

    const auto& vol_field = field()->mesh()->cellsVolumeVector();
    const f64 sign = (SourceSign == Sign::Positive) ? -1.0 : 1.0;
    for (Eigen::Index i = 0; i < vol_field.size(); ++i) {
        insert(static_cast<size_t>(i), static_cast<size_t>(i), sign * _coeff * vol_field[i]);
    }
}
} // namespace prism::scheme::source
