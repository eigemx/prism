#include "model.h"
#include "prism/equation/transport.h"
#include "prism/field/density.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::turbulence {

class KOmega : public IRASModel {
  public:
    KOmega(const SharedPtr<field::Density>& rho,
           const SharedPtr<field::Velocity>& mdot,
           const SharedPtr<field::Velocity>& U,
           const SharedPtr<field::Scalar>& mu,
           f64 dt = 0.0);

    auto k() -> const SharedPtr<field::Scalar>&;
    auto omega() -> const SharedPtr<field::Scalar>&;
    auto turbulentViscosity() -> const SharedPtr<field::Scalar>&;

    auto kEqn() -> eqn::Transport;
    auto omegaEqn() -> eqn::Transport;

  private:
    auto kProduction() -> SharedPtr<field::Scalar>;
    void setViscosity();

    SharedPtr<field::Scalar> _k, _omega, _mu, _mu_t, _mu_eff;
    SharedPtr<field::Velocity> _U, _mdot;
    SharedPtr<field::Density> _rho;
    f64 _c_mu = 0.09;
    f64 _c_alpha1 = 5.0 / 9.0;
    f64 _c_beta1 = 0.075;
    f64 _beta_star = 0.09;
    f64 _sigma_k1 = 2;
    f64 _sigma_w1 = 2;
    f64 _pr_t = 0.9;
    f64 _dt = 0.0;
    bool _is_transient = false;
};
} // namespace prism::turbulence
