#pragma once

#include "prism/field/scalar.h"
#include "prism/operations/operations.h"
#include "source.h"

namespace prism::scheme::source {
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

template <Sign SourceSign, typename Kappa, typename Field>
Laplacian<SourceSign, Kappa, Field>::Laplacian(Kappa kappa, Field phi)
    : IExplicitSource(phi.mesh()->cellCount()), _kappa(kappa), _phi(phi) {}

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


} // namespace prism::scheme::source
