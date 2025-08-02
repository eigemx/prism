#pragma once

#include <optional>
#include <stdexcept>

#include "prism/field/density.h"
#include "prism/field/ifield.h"
#include "prism/log.h"
#include "prism/scheme/temporal/backward_euler.h"
#include "temporal.h"

namespace prism::scheme::temporal {

template <field::IScalarBased Field>
class AdamMoulton : public ITemporal<Field> {
  public:
    AdamMoulton(Field phi);
    AdamMoulton(Field phi, double dt);
    AdamMoulton(field::Density rho, Field phi);
    AdamMoulton(field::Density rho, Field phi, double dt);

    void apply() override;
    auto field() -> Field override;

    auto timeStep() const noexcept -> double;
    void setTimeStep(double dt) noexcept;

    using FieldType = Field;

  private:
    void applyIncompressible();
    void applyCompressible();
    void applyBackwardEuler();

    Field _phi;
    Optional<field::Density> _rho = std::nullopt;
    double _dt {1e-5}; // To avoid initialization with 0
    std::size_t _n_timesteps {0};
};

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(Field phi)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi) initialized for field `{}` with default dt = "
        "{}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(Field phi, double dt)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _dt(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::AdamMoulton(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi, dt) initialized for field `{}` with dt = {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(field::Density rho, Field phi)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(field::Density rho, Field phi, double dt)
    : ITemporal<Field>(phi.mesh()->cellCount()), _phi(phi), _rho(rho), _dt(dt) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        phi.name(),
        _dt);
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::apply() {
    if (_rho.has_value()) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
    _n_timesteps++;
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyCompressible() {
    throw std::runtime_error(
        "prism::scheme::temporal::AdamMoulton::applyCompressible() not implemented yet.");
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyIncompressible() {
    log::debug("Applying Adam-Moulton scheme to field `{}` (no density field provided)",
               _phi.name());

    if (_n_timesteps == 0) {
        // fall back to backward Euler first order scheme, because we have only one time step in
        // the past, and Adam-Moulton requires at least two time steps.
        log::debug(
            "prism::scheme::temporal::AdamMoulton::applyIncompressible() falling back to "
            "Backward-Euler transient scheme for the first time step");
        applyBackwardEuler();
        return;
    }
    const auto& vol_field = _phi.mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!_phi.prevValues().has_value() || !_phi.prevPrevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::AdamMoulton::apply() was called for field `{}`, but "
            "the field does not have previous time steps values stored. Adam-Moulton transient "
            "scheme requires two previous time steps. Please make sure that you set "
            "history tracking of the field using field.setHistorySize(steps) to at least 2.",
            _phi.name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = _phi.prevValues().value();
    const VectorXd phi_prev_prev = _phi.prevPrevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    this->matrix().setIdentity();
    this->matrix().diagonal() = (3.0 / 2.0) * vol_field / _dt;
    this->rhs() = 2.0 * vol_field.cwiseProduct(phi_prev) / _dt;
    this->rhs() += -0.5 * vol_field.cwiseProduct(phi_prev_prev) / _dt;
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyBackwardEuler() {
    BackwardEuler<Field> b_euler(_phi, _dt);
    b_euler.apply();
    this->matrix() = b_euler.matrix();
    this->rhs() = b_euler.rhs();
}

template <field::IScalarBased Field>
auto AdamMoulton<Field>::timeStep() const noexcept -> double {
    return _dt;
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::setTimeStep(double dt) noexcept {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::AdamMoulton::setTimeStep(): dt must be positive");
    }
    _dt = dt;
}

template <field::IScalarBased Field>
auto AdamMoulton<Field>::field() -> Field {
    return _phi;
}

} // namespace prism::scheme::temporal
