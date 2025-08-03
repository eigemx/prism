#pragma once

#include "prism/field/ifield.h"
#include "source.h"

namespace prism::scheme::source {
template <Sign SourceSign, field::IScalarBased Field>
class ImplicitField : public IFullScheme<Field>, public IImplicitSource {
  public:
    ImplicitField(Field& phi) : IFullScheme<Field>(phi) {}
    ImplicitField(double coeff, Field& phi) : IFullScheme<Field>(phi), _coeff(coeff) {}

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    void applyInterior(const mesh::Face& face) override {}
    void applyBoundary() override {}

    double _coeff {1.0};
};


template <Sign SourceSign, field::IScalarBased Field>
void ImplicitField<SourceSign, Field>::apply() {
    this->matrix().setIdentity();
    this->matrix().diagonal() *= _coeff * this->field().mesh()->cellsVolumeVector();

    if constexpr (SourceSign == Sign::Positive) {
        this->matrix() *= -1;
        return;
    }
}
} // namespace prism::scheme::source
