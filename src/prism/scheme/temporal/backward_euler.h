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

    using FieldType = Field;

  private:
    void applyIncompressible();
    void applyCompressible();

    Optional<field::Density> _rho = std::nullopt;
};

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(Field phi) : ITemporal<Field>(phi) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler() initialized for field `{}` with default dt = "
        "{}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(Field phi, double dt) : ITemporal<Field>(phi, dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::BackwardEuler(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::BackwardEuler(phi, dt) initialized for field `{}` with dt = {}",
        this->field().name(),
        this->timeStep());
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi)
    : ITemporal<Field>(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        this->field().name(),
        this->timestep());
}

template <field::IScalarBased Field>
BackwardEuler<Field>::BackwardEuler(field::Density rho, Field phi, double dt)
    : ITemporal<Field>(phi, dt), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        this->field().name(),
        this->timestep());
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::apply() {
    if (_rho.has_value()) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
    this->incrementTimestep();
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::applyCompressible() {
    throw std::runtime_error(
        "prism::scheme::temporal::BackwardEuler::applyCompressible() not implemented yet.");
}

template <field::IScalarBased Field>
void BackwardEuler<Field>::applyIncompressible() {
    log::debug("Applying backward Euler scheme to field `{}` (no density field provided)",
               this->field().name());
    const auto& vol_field = this->field().mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!this->field().prevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::BackwardEuler::apply() was called for field `{}`, but "
            "the field does not have any previous time step values stored.",
            this->field().name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = this->field().prevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    this->matrix().setIdentity();
    this->matrix().diagonal() = vol_field / this->timeStep();
    this->rhs() = vol_field.cwiseProduct(phi_prev) / this->timeStep();
}

} // namespace prism::scheme::temporal
