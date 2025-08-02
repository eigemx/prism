#pragma once

#include <concepts>

#include "prism/scheme/scheme.h"

namespace prism::scheme::temporal {

template <typename Field>
class ITemporal : public scheme::IFullScheme<Field> {
  public:
    ITemporal(std::size_t n_cells) : scheme::IFullScheme<Field>(n_cells) {}

    auto needsCorrection() const noexcept -> bool override;

  private:
    void applyInterior(const mesh::Face& face) override {}
    void applyBoundary() override {}
};


template <typename T>
concept ITemporalBased = std::derived_from<T, ITemporal<typename T::FieldType>>;


template <typename Field>
auto ITemporal<Field>::needsCorrection() const noexcept -> bool {
    return true;
}

} // namespace prism::scheme::temporal
