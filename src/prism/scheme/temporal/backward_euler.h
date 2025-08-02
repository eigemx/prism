#pragma once

#include <optional>
#include <stdexcept>

#include "prism/field/density.h"
#include "prism/field/ifield.h"
#include "prism/log.h"
#include "temporal.h"

namespace prism::scheme::temporal {

template <field::IScalarBased Field>
class BackwardEuler : public ITemporal<Field> {
  public:
    BackwardEuler(Field phi);
    BackwardEuler(Field phi, double dt);
    BackwardEuler(field::Density rho, Field phi);
    BackwardEuler(field::Density rho, Field phi, double dt);

    void apply() override;
    auto field() -> Field override;

    auto timeStep() const noexcept -> double;
    void setTimeStep(double dt) noexcept;

    using FieldType = Field;

  private:
    void applyIncompressible();
    void applyCompressible();

    Field _phi;
    Optional<field::Density> _rho = std::nullopt;
    double _dt {1e-5}; // To avoid initialization with 0
    std::size_t _n_timesteps {0};
};

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(Field phi)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(phi) initialized for field `{}` with default dt = "
        "{}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(Field phi, double dt)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _dt(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::BackwardEuler(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::BackwardEuler(phi, dt) initialized for field `{}` with dt = {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi, double dt)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _rho(rho), _dt(dt) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::apply() {
    if (_rho.has_value()) {
        applyCompressible();
        return;
    }
    applyIncompressible();
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::applyCompressible() {
    throw std::runtime_error(
        "prism::scheme::temporal::BackwardEuler::applyCompressible() not implemented yet.");
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::applyIncompressible() {
    log::debug("Applying backward Euler scheme to field `{}` (no density field provided)",
               _phi.name());
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!_phi.prevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::BackwardEuler::apply() was called for field `{}`, but "
            "the field does not have any previous time step values stored.",
            _phi.name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = _phi.prevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    this->matrix().setIdentity();
    this->matrix().diagonal() = vol_field / _dt;
    this->rhs() = vol_field.cwiseProduct(phi_prev) / _dt;
}

template <field::IScalarBased Field>
auto BackwardEuler<Field>::timeStep() const noexcept -> double {
    return _dt;
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::setTimeStep(double dt) noexcept {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::BackwardEuler::setTimeStep(): dt must be positive");
    }
    _dt = dt;
}

template <field::IScalarBased Field>
auto BackwardEuler<Field>::field() -> Field {
    return _phi;
}

} // namespace prism::scheme::temporal
