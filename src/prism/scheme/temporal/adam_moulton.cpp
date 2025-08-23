#include "adam_moulton.h"

#include "prism/log.h"
#include "prism/scheme/temporal/backward_euler.h"
#include <stdexcept>

namespace prism::scheme::temporal {

AdamMoulton::AdamMoulton(const SharedPtr<field::Scalar>& phi) : ITemporal(phi) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi) initialized for field `{}` with default dt = "
        "{}",
        field()->name(),
        timeStep());
}

AdamMoulton::AdamMoulton(const SharedPtr<field::Scalar>& phi, double dt) : ITemporal(phi, dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::AdamMoulton(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::AdamMoulton(phi, dt) initialized for field `{}` with dt = {}",
        field()->name(),
        timeStep());
}

AdamMoulton::AdamMoulton(const SharedPtr<field::Density>& rho,
                         const SharedPtr<field::Scalar>& phi)
    : ITemporal(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        field()->name(),
        timeStep());
}

AdamMoulton::AdamMoulton(const SharedPtr<field::Density>& rho,
                         const SharedPtr<field::Scalar>& phi,
                         double dt)
    : ITemporal(phi, dt), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::AdamMoulton(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        field()->name(),
        timeStep());
}

void AdamMoulton::apply() {
    if (_rho) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
    incrementTimestep();
}

void AdamMoulton::applyCompressible() { // NOLINT
    throw std::runtime_error(
        "prism::scheme::temporal::AdamMoulton::applyCompressible() not implemented yet.");
}

void AdamMoulton::applyIncompressible() {
    log::debug("Applying Adam-Moulton scheme to field `{}` (no density field provided)",
               field()->name());

    if (timeStepsCount() == 0) {
        // fall back to backward Euler first order scheme, because we have only one time step in
        // the past, and Adam-Moulton requires at least two time steps.
        log::debug(
            "prism::scheme::temporal::AdamMoulton::applyIncompressible() falling back to "
            "Backward-Euler transient scheme for the first time step");
        applyBackwardEuler();
        return;
    }
    const auto& vol_field = field()->mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!field()->prevValues().has_value() || !field()->prevPrevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::AdamMoulton::apply() was called for field `{}`, but "
            "the field does not have previous time steps values stored. Adam-Moulton transient "
            "scheme requires two previous time steps. Please make sure that you set "
            "history tracking of the field using field.setHistorySize(steps) to at least 2.",
            field()->name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = field()->prevValues().value();
    const VectorXd phi_prev_prev = field()->prevPrevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    matrix().setIdentity();
    matrix().diagonal() = (3.0 / 2.0) * vol_field / timeStep();
    rhs() = 2.0 * vol_field.cwiseProduct(phi_prev) / timeStep();
    rhs() += -0.5 * vol_field.cwiseProduct(phi_prev_prev) / timeStep();
}

void AdamMoulton::applyBackwardEuler() {
    BackwardEuler b_euler(field(), timeStep());
    b_euler.apply();
    matrix() = b_euler.matrix();
    rhs() = b_euler.rhs();
}

} // namespace prism::scheme::temporal
