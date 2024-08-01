#pragma once

#include <spdlog/spdlog.h>

#include "prism/constants.h"
#include "prism/equation/equation.h"

namespace prism::solver {

template <typename Field>
class RelaxerBase {
  public:
    // pre_relax() is called before solving the equation (suitable for implicit schemes)
    virtual void pre_relax(TransportEquation<Field>& eqn, double lambda) const = 0;

    // post_relax() is called after solving the equation (suitable for explicit schemes)
    virtual void post_relax(TransportEquation<Field>& eqn, double lambda) const = 0;
};

template <typename Field>
class ImplicitUnderRelaxation : public RelaxerBase<Field> {
  public:
    void pre_relax(TransportEquation<Field>& eqn, double lambda) const override;

    // no post_relax() is needed for implicit under-relaxation
    void inline post_relax(TransportEquation<Field>& eqn, double lambda) const override {}
};

template <typename Field>
class ExplicitUnderRelaxation : public RelaxerBase<Field> {
  public:
    // no pre_relax() is needed for explicit under-relaxation
    void inline pre_relax(TransportEquation<Field>& eqn, double lambda) const override {}

    void post_relax(TransportEquation<Field>& eqn, double lambda) const override;
};

template <typename Field>
void ImplicitUnderRelaxation<Field>::pre_relax(TransportEquation<Field>& eqn,
                                               double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        spdlog::debug(
            "ImplicitUnderRelaxation::pre_relax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& A = eqn.matrix();
    auto& b = eqn.rhs();
    const auto& phi = eqn.field().data();

    // Moukallad et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    A.diagonal() /= lambda;
    b += (1 - lambda) * A.diagonal().cwiseProduct(phi);
}

template <typename Field>
void ExplicitUnderRelaxation<Field>::post_relax(TransportEquation<Field>& eqn,
                                                double lambda) const {
    if (lambda < 0.0 || lambda == 1.0) {
        spdlog::debug(
            "ImplicitUnderRelaxation::pre_relax(): relaxation factor is"
            " either less than 0.0 or equals 1.0, skipping relaxation...");
        return;
    }

    auto& phi = eqn.field().data();
    const auto& phi_old = eqn.field_prev_iter().data();

    // Moukallad et. al, 14.2.2  Explicit Under-Relaxation
    // equation 14.6
    phi = phi_old + (lambda * (phi - phi_old));
}

} // namespace prism::solver