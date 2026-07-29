#include "KOmega.h"

#include "prism/scheme/source/constant_scalar.h"
#include "prism/scheme/temporal/adam_moulton.h"

namespace prism::turbulence {
KOmega::KOmega(const SharedPtr<field::Density>& rho,
               const SharedPtr<field::Velocity>& mdot,
               const SharedPtr<field::Velocity>& U,
               const SharedPtr<field::Scalar>& mu,
               f64 dt)
    : _rho(rho), _mdot(mdot), _U(U), _mu(mu), _dt(dt) {
    _k = makeShared<field::Scalar>("k", mu->mesh(), 0.0);
    _omega = makeShared<field::Scalar>("omega", mu->mesh(), 0.0);
    _mu_t = makeShared<field::Scalar>("mu_t", mu->mesh(), 0.0);
    _mu_eff = makeShared<field::Scalar>("mu_eff", mu->mesh(), 0.0);

    if (dt > 0.0) {
        _is_transient = true;
    }
}

void KOmega::setViscosity() {
    auto mu_t_vals = field::mul(_rho, field::divide(_k, _omega));
    _mu_t->values() = mu_t_vals;
    _mu_eff->values() = mu_t_vals + _mu->values();
}

auto KOmega::equations() -> Pair<eqn::Transport, eqn::Transport> {
    setViscosity();
    auto pk = kProduction();
    return {buildKEqn(pk), buildOmegaEqn(pk)};
}

auto KOmega::buildKEqn(const SharedPtr<field::Scalar>& pk) -> eqn::Transport {
    auto sk_values = field::mul(_rho, _k, _omega);
    auto sk = makeShared<field::Scalar>("sk", _k->mesh(), _beta_star * sk_values);
    auto k_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _k),
                                scheme::diffusion::NonCorrected(_mu_eff, _k),
                                scheme::source::ConstantScalar<Sign::Positive>(pk),
                                scheme::source::ConstantScalar<Sign::Negative>(sk));

    if (_is_transient) {
        k_eqn.addScheme(scheme::temporal::AdamMoulton(_rho, _k, _dt));
    }

    return k_eqn;
}

auto KOmega::buildOmegaEqn(const SharedPtr<field::Scalar>& pk) -> eqn::Transport {
    VectorXd sk_values = _c_alpha1 * field::mul(field::divide(_omega, _k), pk);
    sk_values -= _c_beta1 * field::mul(_rho, _omega, _omega);
    auto sk = makeShared<field::Scalar>("sk", _k->mesh(), _beta_star * sk_values);
    auto omega_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _omega),
                                    scheme::diffusion::NonCorrected(_mu_eff, _omega),
                                    scheme::source::ConstantScalar<Sign::Positive>(sk));

    if (_is_transient) {
        omega_eqn.addScheme(scheme::temporal::AdamMoulton(_rho, _omega, _dt));
    }

    return omega_eqn;
}
} // namespace prism::turbulence
