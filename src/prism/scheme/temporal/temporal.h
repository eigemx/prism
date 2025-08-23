#pragma once

#include <concepts>

#include "prism/scheme/scheme.h"

namespace prism::scheme::temporal {

class ITemporal : public scheme::IFullScheme {
  public:
    ITemporal(const SharedPtr<field::Scalar>& phi);
    ITemporal(const SharedPtr<field::Scalar>& phi, double dt);

    auto needsCorrection() const noexcept -> bool override;

    auto timeStep() const noexcept -> double;
    void setTimeStep(double dt);

  protected:
    auto timeStepsCount() const noexcept -> std::size_t;
    void incrementTimestep() noexcept;
    void resetTimestep() noexcept;

  private:
    void applyInterior(const mesh::Face& face) override {}
    void applyBoundary() override {}

    double _dt {1e-5}; // To avoid initialization with 0
    std::size_t _n_timesteps {0};
};


template <typename T>
concept ITemporalBased = std::derived_from<T, ITemporal>;

} // namespace prism::scheme::temporal
