#include "KEpsilon.h"

#include <memory>

#include "k_production.h"
#include "prism/scheme/source/constant_scalar.h"
#include "prism/scheme/temporal/adam_moulton.h"

namespace prism::turbulence {

IncompressibleKEpsilon::IncompressibleKEpsilon(const SharedPtr<field::Velocity>& mdot,
                                               const SharedPtr<field::Velocity>& U,
                                               const SharedPtr<field::Scalar>& mu,
                                               Optional<f64> dt)
    : _mdot(mdot), _U(U), _nu(mu) {
    _k = std::make_shared<field::Scalar>("k", mu->mesh(), 0.0);
    _epsilon = std::make_shared<field::Scalar>("epsilon", mu->mesh(), 0.0);
    _nu_t = std::make_shared<field::Scalar>("mu_t", mu->mesh(), 0.0);

    if (dt.has_value() && dt.value() > 0.0) {
        _is_transient = true;
        _dt = dt.value();
        _k->setHistorySize(2);
        _epsilon->setHistorySize(2);
    }
}

void IncompressibleKEpsilon::setTurbulentViscosity() {
    VectorXd _nu_t_vals = _c_mu * field::divide(field::mul(_k, _k), _epsilon);
    _nu_t->values() = _nu_t_vals;
}

auto IncompressibleKEpsilon::kEqn() -> eqn::Transport {
    setTurbulentViscosity();
    auto pk = kProductionIncompressible(_nu_t, _k, _U);

    VectorXd nu_eff_k_vals = _nu->values() + (_nu_t->values() / _sigma_k);
    auto nu_eff_k = std::make_shared<field::Scalar>("mu_eff_k", _k->mesh(), nu_eff_k_vals);

    auto k_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _k),
                                scheme::diffusion::NonCorrected(nu_eff_k, _k),
                                scheme::source::ConstantScalar<Sign::Positive>(pk),
                                scheme::source::ConstantScalar<Sign::Negative>(_epsilon));

    if (_is_transient) {
        k_eqn.addScheme(scheme::temporal::AdamMoulton(_k, _dt));
    }

    return k_eqn;
}

auto IncompressibleKEpsilon::epsilonEqn() -> eqn::Transport {
    setTurbulentViscosity();
    auto pk = kProductionIncompressible(_nu_t, _k, _U);

    VectorXd nu_eff_eps_vals = _nu->values() + (_nu_t->values() / _sigma_eps);
    auto nu_eff_eps = std::make_shared<field::Scalar>("mu_eff_eps", _k->mesh(), nu_eff_eps_vals);

    auto eps_over_k = field::divide(_epsilon, _k);

    auto eps_source = std::make_shared<field::Scalar>(
        "eps_source", _k->mesh(), _c_eps1 * field::mul(eps_over_k, pk));

    auto eps_sink = std::make_shared<field::Scalar>(
        "eps_sink", _k->mesh(), _c_eps2 * field::mul(_epsilon, eps_over_k));

    auto eps_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _epsilon),
                                  scheme::diffusion::NonCorrected(nu_eff_eps, _epsilon),
                                  scheme::source::ConstantScalar<Sign::Positive>(eps_source),
                                  scheme::source::ConstantScalar<Sign::Negative>(eps_sink));

    if (_is_transient) {
        eps_eqn.addScheme(scheme::temporal::AdamMoulton(_epsilon, _dt));
    }

    return eps_eqn;
}

auto IncompressibleKEpsilon::k() -> const SharedPtr<field::Scalar>& {
    return _k;
}

auto IncompressibleKEpsilon::epsilon() -> const SharedPtr<field::Scalar>& {
    return _epsilon;
}

auto IncompressibleKEpsilon::turbulentViscosity() -> const SharedPtr<field::Scalar>& {
    return _nu_t;
}

} // namespace prism::turbulence
