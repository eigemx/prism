#pragma once

#include "prism/equation/transport.h"
#include "prism/log.h"

namespace prism::solver {

template <typename Field>
class IRelaxer {
  public:
    // preRelax() is called before solving the equation (suitable for implicit schemes)
    virtual void preRelax(eqn::Transport<Field>& eqn, double lambda) const = 0;

    // postRelax() is called after solving the equation (suitable for explicit schemes)
    virtual void postRelax(eqn::Transport<Field>& eqn, double lambda) const = 0;
};

template <typename Field>
class ImplicitUnderRelaxation : public IRelaxer<Field> {
  public:
    void preRelax(eqn::Transport<Field>& eqn, double lambda) const override;

    // no post_relax() is needed for implicit under-relaxation
    void inline postRelax(eqn::Transport<Field>& eqn, double lambda) const override {}
};

template <typename Field>
class ExplicitUnderRelaxation : public IRelaxer<Field> {
  public:
    // no pre_relax() is needed for explicit under-relaxation
    void inline preRelax(eqn::Transport<Field>& eqn, double lambda) const override {}

    void postRelax(eqn::Transport<Field>& eqn, double lambda) const override;
};

template <typename Field>
void ImplicitUnderRelaxation<Field>::preRelax(eqn::Transport<Field>& eqn, double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        log::debug(
            "ImplicitUnderRelaxation::preRelax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& A = eqn.matrix();
    auto& b = eqn.rhs();
    const auto& phi = eqn.field().values();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= lambda;
    b += (1 - lambda) * A.diagonal().cwiseProduct(phi);
}

template <typename Field>
void ExplicitUnderRelaxation<Field>::postRelax(eqn::Transport<Field>& eqn, double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        log::debug(
            "ImplicitUnderRelaxation::preRelax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& phi = eqn.field().values();
    const auto& phi_old = eqn.field_prev_iter().values();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda * (phi - phi_old));
}

} // namespace prism::solver