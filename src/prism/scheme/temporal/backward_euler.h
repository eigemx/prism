#pragma once

#include <optional>
#include <stdexcept>

#include "prism/field/density.h"
#include "prism/log.h"
#include "temporal.h"

namespace prism::scheme::temporal {

template <typename Field>
class BackwardEuler : public ITemporal<Field> {
  public:
    BackwardEuler(Field phi);
    BackwardEuler(Field phi, double dt);
    BackwardEuler(field::Density rho, Field phi);
    BackwardEuler(field::Density rho, Field phi, double dt);

    void apply() override;
    auto field() -> Field override;

    auto timeStep() const -> double;
    void setTimeStep(double dt);

    using FieldType = Field;

  private:
    Field _phi;
    std::optional<field::Density> _rho = std::nullopt;
    double _dt {1e-5}; // To avoid initialization with 0
};

template <typename Field>
BackwardEuler<Field>::BackwardEuler(Field phi) : _phi(phi) {
    log::debug(
        "scheme::temporal::BackwardEuler() initialized for field `{}` with default dt = {}",
        phi.name(),
        _dt);
}

template <typename Field>
BackwardEuler<Field>::BackwardEuler(Field phi, double dt) : _phi(phi), _dt(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("scheme::temporal::BackwardEuler(): dt must be positive");
    }
    log::debug("scheme::temporal::BackwardEuler() initialized for field `{}` with dt = {}",
               phi.name(),
               _dt);
}

template <typename Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi)
    : BackwardEuler(rho, phi, 1e-5), _rho(rho) {
    log::debug(
        "scheme::temporal::BackwardEuler() initialized for field `{}` with default dt = {}",
        phi.name(),
        _dt);
}

template <typename Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi, double dt)
    : BackwardEuler(phi, dt), _rho(rho) {
    log::debug("scheme::temporal::BackwardEuler() initialized for field `{}` with dt = {}",
               phi.name(),
               _dt);
}

template <typename Field>
void BackwardEuler<Field>::apply() {
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();
    this->matrix().setIdentity();
    this->matrix().diagonal() = vol_field / _dt;
    this->rhs() = - vol_field.cwiseProduct(_phi.values()) / _dt;
}

template <typename Field>
auto BackwardEuler<Field>::timeStep() const -> double {
    return _dt;
}

template <typename Field>
void BackwardEuler<Field>::setTimeStep(double dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("scheme::temporal::BackwardEuler(): dt must be positive");
    }
    _dt = dt;
}
} // namespace prism::scheme::temporal
