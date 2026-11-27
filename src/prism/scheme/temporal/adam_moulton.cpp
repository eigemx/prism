#include "adam_moulton.h"

#include <stdexcept>

#include "prism/log.h"
#include "prism/scheme/temporal/backward_euler.h"

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

AdamMoulton::AdamMoulton(const SharedPtr<field::Density>& rho, const SharedPtr<field::Scalar>& phi)
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
    resetValues();
    if (_rho) {
        applyCompressible();
    } else {
        applyIncompressible();
    }
}

void AdamMoulton::applyCompressible() { // NOLINT
    throw std::runtime_error(
        "prism::scheme::temporal::AdamMoulton::applyCompressible() not implemented yet.");
}

void AdamMoulton::applyIncompressible() {
    log::debug("Applying Adam-Moulton scheme to field `{}` (no density field provided)",
               field()->name());

    // we need to make sure that the field is keeping track of at least one time step
    auto phi_prev_opt = field()->prevValues();
    if (!phi_prev_opt) {
        throw std::runtime_error(fmt::format(
            "prism::scheme::temporal::AdamMoulton::apply() was called for field `{}`, but "
            "the field does not have previous time steps values stored. Adam-Moulton transient "
            "scheme requires two previous time steps. Please make sure that you set "
            "history tracking of the field using field.setHistorySize(steps) to at least 2.",
            field()->name()));
    }

    auto phi_prev_prev_opt = field()->prevPrevValues();
    if (!phi_prev_prev_opt) {
        // fall back to backward Euler first order scheme, because we have only one time step in
        // the past, and Adam-Moulton requires at least two time steps.
        log::debug(
            "prism::scheme::temporal::AdamMoulton::applyIncompressible() falling back to "
            "Backward-Euler transient scheme for the first time step");
        const VectorXd& phi_prev = *phi_prev_opt;
        detail::assembleBackwardEuler(*this, field(), timeStep(), phi_prev);
        return;
    }
    const auto& vol_field = field()->mesh()->cellsVolumeVector();

    const VectorXd& phi_prev = *phi_prev_opt;
    const VectorXd& phi_prev_prev = *phi_prev_prev_opt;

    // Note that the left hand side is constant for all time steps, we need to utilize this to
    // avoid recalculation of the LHS matrix each time step.
    const f64 coeff = 1.5 / timeStep();
    for (Eigen::Index i = 0; i < vol_field.size(); ++i) {
        insert(static_cast<size_t>(i), static_cast<size_t>(i), coeff * vol_field[i]);
    }
    rhs() = 2.0 * vol_field.cwiseProduct(phi_prev) / timeStep();
    rhs() += -0.5 * vol_field.cwiseProduct(phi_prev_prev) / timeStep();
}

} // namespace prism::scheme::temporal
