#include "backward_euler.h"

#include <stdexcept>

#include "prism/log.h"

namespace prism::scheme::temporal {

BackwardEuler::BackwardEuler(const SharedPtr<field::Scalar>& phi) : ITemporal(phi) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler() initialized for field `{}` with default dt = "
        "{}",
        field()->name(),
        timeStep());
}

BackwardEuler::BackwardEuler(const SharedPtr<field::Scalar>& phi, double dt)
    : ITemporal(phi, dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument(
            "prism::scheme::temporal::BackwardEuler(phi, dt): dt must be positive");
    }
    log::debug(
        "prism::scheme::temporal::BackwardEuler(phi, dt) initialized for field `{}` with dt = {}",
        field()->name(),
        timeStep());
}

BackwardEuler::BackwardEuler(const SharedPtr<field::Density>& rho,
                             const SharedPtr<field::Scalar>& phi)
    : ITemporal(phi), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi) initialized for field `{}` with "
        "default dt = {}",
        field()->name(),
        timeStep());
}

BackwardEuler::BackwardEuler(const SharedPtr<field::Density>& rho,
                             const SharedPtr<field::Scalar>& phi,
                             double dt)
    : ITemporal(phi, dt), _rho(rho) {
    log::debug(
        "prism::scheme::temporal::BackwardEuler(rho, phi, dt) initialized for field `{}` with dt "
        "= {}",
        field()->name(),
        timeStep());
}

void BackwardEuler::apply() {
    if (_rho) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
    incrementTimestep();
}

void BackwardEuler::applyCompressible() { // NOLINT
    throw std::runtime_error(
        "prism::scheme::temporal::BackwardEuler::applyCompressible() not implemented yet.");
}

void BackwardEuler::applyIncompressible() {
    log::debug("Applying backward Euler scheme to field `{}` (no density field provided)",
               field()->name());
    const auto& vol_field = field()->mesh()->cellsVolumeVector();

    // we need to make sure that the field is keeping track of at least one time step
    if (!field()->prevValues().has_value()) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::BackwardEuler::apply() was called for field `{}`, but "
            "the field does not have any previous time step values stored.",
            field()->name()));
    }

    /// TODO: when we find a way to avoid the copy, we need to adjust the code below to const&
    const VectorXd phi_prev = field()->prevValues().value();

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    matrix().setIdentity();
    matrix().diagonal() = vol_field / timeStep();
    rhs() = vol_field.cwiseProduct(phi_prev) / timeStep();
}


} // namespace prism::scheme::temporal
