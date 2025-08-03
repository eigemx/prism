#pragma once

#include <concepts>

#include "prism/scheme/scheme.h"

namespace prism::scheme::temporal {

template <field::IScalarBased Field>
class ITemporal : public scheme::IFullScheme<Field> {
  public:
    ITemporal(Field phi);
    ITemporal(Field phi, double dt);

    auto needsCorrection() const noexcept -> bool override;

    auto timeStep() const noexcept -> double;
    void setTimeStep(double dt) noexcept;

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
concept ITemporalBased = std::derived_from<T, ITemporal<typename T::FieldType>>;

template <field::IScalarBased Field>
ITemporal<Field>::ITemporal(Field phi) : scheme::IFullScheme<Field>(phi) {}

template <field::IScalarBased Field>
ITemporal<Field>::ITemporal(Field phi, double dt) : scheme::IFullScheme<Field>(phi), _dt(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::ITemporal(phi, dt): dt must be positive");
    }
}

template <field::IScalarBased Field>
auto ITemporal<Field>::needsCorrection() const noexcept -> bool {
    return true;
}

template <field::IScalarBased Field>
auto ITemporal<Field>::timeStep() const noexcept -> double {
    return _dt;
}

template <field::IScalarBased Field>
void ITemporal<Field>::setTimeStep(double dt) noexcept {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::ITemporal::setTimeStep(): dt must be positive");
    }
    _dt = dt;
}

template <field::IScalarBased Field>
auto ITemporal<Field>::timeStepsCount() const noexcept -> std::size_t {
    return _n_timesteps;
}

template <field::IScalarBased Field>
void ITemporal<Field>::incrementTimestep() noexcept {
    _n_timesteps++;
}

template <field::IScalarBased Field>
void ITemporal<Field>::resetTimestep() noexcept {
    _n_timesteps = 0;
}

} // namespace prism::scheme::temporal
