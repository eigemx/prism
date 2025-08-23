#include "transport.h"

namespace prism::eqn {
void Transport::updateCoeffs() {
    // iterate over all equation's finite volume schemes
    for (auto& scheme : _schemes) {
        // apply the scheme
        scheme->apply();

        // update equation's universal coefficient matrix and RHS vector
        matrix() += scheme->matrix();
        rhs() += scheme->rhs();
    }

    for (auto& source : _sources) {
        source->apply();
        rhs() += source->rhs();
    }

    prism::boundary::detail::applyBoundaryIfExists("prism::eqn::Transport", *this);
}

void Transport::zeroOutCoeffs() {
    for (auto& scheme : _schemes) {
        scheme->matrix().setZero();
        scheme->rhs().setZero();
    }

    for (auto& source : _sources) {
        source->rhs().setZero();
    }

    // zero out the universal coefficient matrix and RHS vector
    matrix().setZero();
    rhs().setZero();
}

auto Transport::underRelaxFactor() const -> double {
    return _relaxation_factor;
}

void Transport::setUnderRelaxFactor(double factor) {
    if (factor < 0.0 || factor > 1.0) {
        log::warn(
            "Transport::setUnderRelaxFactor(): relaxation factor must be in [0, 1] range. "
            "Ignoring the value {} and using 1.0 instead.",
            factor);
        factor = 1.0;
        return;
    }
    _relaxation_factor = factor;
}

void Transport::relax() {
    auto& A = matrix();
    auto& b = rhs();
    const auto& phi = field()->values();

    // Moukalled et. al, 14.2 Under-Relaxation of the Algebraic Equations, equation 14.9.
    log::debug("Transport::relax(): applying implicit under-relaxation with factor = {}",
               _relaxation_factor);
    b += ((1.0 - _relaxation_factor) / _relaxation_factor) * A.diagonal().cwiseProduct(phi);
    A.diagonal() /= _relaxation_factor;
}

auto Transport::field() const -> const SharedPtr<field::Scalar>& {
    return _phi;
}

auto Transport::convectionScheme() -> SharedPtr<scheme::IFullScheme> {
    return _conv_scheme;
}

auto Transport::diffusionScheme() -> SharedPtr<scheme::IFullScheme> {
    return _diff_scheme;
}

auto Transport::isTransient() const noexcept -> bool {
    return _is_transient;
}

} // namespace prism::eqn
