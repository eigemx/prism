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
    using FieldType = Field;

  private:
    void applyIncompressible();
    void applyCompressible();
    void applyBackwardEuler();

    Optional<field::Density> _rho = std::nullopt;
};

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(Field phi) : ITemporal<Field>(phi) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi) initialized for field `{}` with default dt = "
        "{}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(Field phi, double dt) : ITemporal<Field>(phi, dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::AdamMoulton(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi, dt) initialized for field `{}` with dt = {}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(field::Density rho, Field phi)
    : ITemporal<Field>(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
AdamMoulton<Field>::AdamMoulton(field::Density rho, Field phi, double dt)
    : ITemporal<Field>(phi, dt), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::apply() {
    if (_rho.has_value()) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
    this->incrementTimestep();
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyCompressible() {
    throw std::runtime_error(
        "prism::scheme::temporal::AdamMoulton::applyCompressible() not implemented yet.");
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyIncompressible() {
    log::debug("Applying Adam-Moulton scheme to field `{}` (no density field provided)",
               this->field().name());

    if (this->timeStepsCount() == 0) {
        // fall back to backward Euler first order scheme, because we have only one time step in
        // the past, and Adam-Moulton requires at least two time steps.
        log::debug(
            "prism::scheme::temporal::AdamMoulton::applyIncompressible() falling back to "
            "Backward-Euler transient scheme for the first time step");
        applyBackwardEuler();
        return;
    }
    const auto& vol_field = this->field().mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!this->field().prevValues().has_value() || !this->field().prevPrevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::AdamMoulton::apply() was called for field `{}`, but "
            "the field does not have previous time steps values stored. Adam-Moulton transient "
            "scheme requires two previous time steps. Please make sure that you set "
            "history tracking of the field using field.setHistorySize(steps) to at least 2.",
            this->field().name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = this->field().prevValues().value();
    const VectorXd phi_prev_prev = this->field().prevPrevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    this->matrix().setIdentity();
    this->matrix().diagonal() = (3.0 / 2.0) * vol_field / this->timeStep();
    this->rhs() = 2.0 * vol_field.cwiseProduct(phi_prev) / this->timeStep();
    this->rhs() += -0.5 * vol_field.cwiseProduct(phi_prev_prev) / this->timeStep();
}

template <field::IScalarBased Field>
void AdamMoulton<Field>::applyBackwardEuler() {
    BackwardEuler<Field> b_euler(this->field(), this->timeStep());
    b_euler.apply();
    this->matrix() = b_euler.matrix();
    this->rhs() = b_euler.rhs();
}

} // namespace prism::scheme::temporal
