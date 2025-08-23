#pragma once

#include "prism/field/scalar.h"
#include "prism/operations/operations.h"
#include "source.h"

namespace prism::scheme::source {
template <Sign SourceSign = Sign::Positive>
class Laplacian : public IExplicitSource {
  public:
    Laplacian(const SharedPtr<field::Scalar>& kappa, const SharedPtr<field::Scalar>& phi);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    SharedPtr<field::Scalar> _kappa;
    SharedPtr<field::Scalar> _phi;
};

template <Sign SourceSign>
Laplacian<SourceSign>::Laplacian(const SharedPtr<field::Scalar>& kappa,
                                 const SharedPtr<field::Scalar>& phi)
    : IExplicitSource(phi->mesh()->cellCount()), _kappa(kappa), _phi(phi) {}

template <Sign SourceSign>
void inline Laplacian<SourceSign>::apply() {
    auto grad_phi = ops::grad(*_phi);
    auto div = ops::div(grad_phi);

    if constexpr (SourceSign == Sign::Positive) {
        rhs() = div.values();
        return;
    }
    rhs() = -div.values();
}


} // namespace prism::scheme::source
