#pragma once

#include "prism/field/scalar.h"
#include "prism/log.h"
#include "prism/operations/operations.h"
#include "source.h"

namespace prism::scheme::source {
// Adds a source for a gradient of a scalar field, in a specific coordinate
// for example the gradient of the pressure in the x-direction: ∂p/∂x
template <Sign SourceSign = Sign::Positive, typename Field = field::Scalar>
class Gradient : public IExplicitSource {
  public:
    Gradient(Field phi, Coord coord);
    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return true; }

  private:
    Field _phi;
    Coord _coord;
};

template <Sign SourceSign, typename Field>
Gradient<SourceSign, Field>::Gradient(Field phi, Coord coord)
    : _phi(phi), _coord(coord), IExplicitSource(phi.mesh()->cellCount()) {
    log::debug(
        "prism::scheme::source::Gradient(): Creating {}-coordinate gradient source for field "
        "'{}'",
        field::coordToStr(coord),
        phi.name());
}

template <Sign SourceSign, typename Field>
void Gradient<SourceSign, Field>::apply() {
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();
    rhs() = ops::grad(_phi, _coord).values().cwiseProduct(vol_field);

    if constexpr (SourceSign == Sign::Negative) {
        rhs() = -rhs();
    }
}


} // namespace prism::scheme::source
