#include "temporal.h"


namespace prism::scheme::temporal {
ITemporal::ITemporal(const SharedPtr<field::Scalar>& phi) : scheme::IFullScheme(phi) {}

ITemporal::ITemporal(const SharedPtr<field::Scalar>& phi, double dt)
    : scheme::IFullScheme(phi), _dt(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::ITemporal(phi, dt): dt must be positive");
    }
}

auto ITemporal::needsCorrection() const noexcept -> bool {
    return true;
}

auto ITemporal::timeStep() const noexcept -> double {
    return _dt;
}

void ITemporal::setTimeStep(double dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::ITemporal::setTimeStep(): dt must be positive");
    }
    _dt = dt;
}

auto ITemporal::timeStepsCount() const noexcept -> std::size_t {
    return _n_timesteps;
}

void ITemporal::incrementTimestep() noexcept {
    _n_timesteps++;
}

void ITemporal::resetTimestep() noexcept {
    _n_timesteps = 0;
}
} // namespace prism::scheme::temporal
