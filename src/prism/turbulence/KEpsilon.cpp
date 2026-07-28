#include "KEpsilon.h"

#include <memory>

#include "prism/operations/operations.h"
#include "prism/scheme/source/constant_scalar.h"
#include "prism/scheme/temporal/adam_moulton.h"

namespace prism::turbulence {

KEpsilon::KEpsilon(const SharedPtr<field::Density>& rho,
                   const SharedPtr<field::Velocity>& mdot,
                   const SharedPtr<field::Velocity>& U,
                   const SharedPtr<field::Scalar>& mu,
                   f64 dt)
    : _rho(rho), _mdot(mdot), _U(U), _mu(mu), _dt(dt) {
    _k = std::make_shared<field::Scalar>("k", mu->mesh(), 0.0);
    _epsilon = std::make_shared<field::Scalar>("epsilon", mu->mesh(), 0.0);
    _mu_t = std::make_shared<field::Scalar>("mu_t", mu->mesh(), 0.0);

    if (dt > 0.0) {
        _is_transient = true;
        _k->setHistorySize(2);
        _epsilon->setHistorySize(2);
    }
}

void KEpsilon::setViscosity() {
    VectorXd mu_t_vals = _c_mu * field::mul(_rho, field::divide(field::mul(_k, _k), _epsilon));
    _mu_t->values() = mu_t_vals;
}

auto KEpsilon::kProduction() -> SharedPtr<field::Scalar> {
    const auto& mesh = _k->mesh();
    VectorXd pk_vals(mesh->cellCount());
    const auto& mu_t_vals = _mu_t->values();

    for (const auto& cell : mesh->cells()) {
        Tensor3d gradU = _U->gradAtCell(cell);
        pk_vals[cell.id()] =
            mu_t_vals[cell.id()] * ops::doubleContraction(ops::dev(ops::twoSymm(gradU)), gradU);
    }

    return std::make_shared<field::Scalar>("Pk", mesh, pk_vals);
}

auto KEpsilon::kEqn() -> eqn::Transport {
    setViscosity();
    auto pk = kProduction();

    VectorXd mu_eff_k_vals = _mu->values() + _mu_t->values() / _sigma_k;
    auto mu_eff_k = std::make_shared<field::Scalar>("mu_eff_k", _k->mesh(), mu_eff_k_vals);

    VectorXd dissipation_vals = field::mul(_rho, _epsilon);
    auto dissipation = std::make_shared<field::Scalar>("dissipation", _k->mesh(), dissipation_vals);

    auto k_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _k),
                                scheme::diffusion::NonCorrected(mu_eff_k, _k),
                                scheme::source::ConstantScalar<Sign::Positive>(pk),
                                scheme::source::ConstantScalar<Sign::Negative>(dissipation));

    if (_is_transient) {
        k_eqn.addScheme(scheme::temporal::AdamMoulton(_rho, _k, _dt));
    }

    return k_eqn;
}

auto KEpsilon::epsilonEqn() -> eqn::Transport {
    setViscosity();
    auto pk = kProduction();

    VectorXd mu_eff_eps_vals = _mu->values() + _mu_t->values() / _sigma_eps;
    auto mu_eff_eps = std::make_shared<field::Scalar>("mu_eff_eps", _k->mesh(), mu_eff_eps_vals);

    auto eps_over_k =
        std::make_shared<field::Scalar>("eps_over_k", _k->mesh(), field::divide(_epsilon, _k));

    auto eps_source = std::make_shared<field::Scalar>(
        "eps_source", _k->mesh(), _c_eps1 * field::mul(eps_over_k, pk));

    auto eps_sink = std::make_shared<field::Scalar>(
        "eps_sink", _k->mesh(), _c_eps2 * field::mul(_rho, _epsilon, eps_over_k));

    auto eps_eqn = eqn::Transport(scheme::convection::LinearUpwind(_mdot, _epsilon),
                                  scheme::diffusion::NonCorrected(mu_eff_eps, _epsilon),
                                  scheme::source::ConstantScalar<Sign::Positive>(eps_source),
                                  scheme::source::ConstantScalar<Sign::Negative>(eps_sink));

    if (_is_transient) {
        eps_eqn.addScheme(scheme::temporal::AdamMoulton(_rho, _epsilon, _dt));
    }

    return eps_eqn;
}

auto KEpsilon::k() -> const SharedPtr<field::Scalar>& {
    return _k;
}

auto KEpsilon::epsilon() -> const SharedPtr<field::Scalar>& {
    return _epsilon;
}

auto KEpsilon::turbulentViscosity() -> const SharedPtr<field::Scalar>& {
    return _mu_t;
}

} // namespace prism::turbulence
