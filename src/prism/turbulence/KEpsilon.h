#pragma once

#include "model.h"
#include "prism/equation/transport.h"
#include "prism/field/density.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::turbulence {

class KEpsilon : public IRASModel {
  public:
    KEpsilon(const SharedPtr<field::Density>& rho,
             const SharedPtr<field::Velocity>& mdot,
             const SharedPtr<field::Velocity>& U,
             const SharedPtr<field::Scalar>& mu,
             f64 dt = 0.0);

    auto k() -> const SharedPtr<field::Scalar>&;
    auto epsilon() -> const SharedPtr<field::Scalar>&;
    auto turbulentViscosity() -> const SharedPtr<field::Scalar>&;

    auto kEqn() -> eqn::Transport;
    auto epsilonEqn() -> eqn::Transport;

  private:
    void setViscosity();
    auto kProduction() -> SharedPtr<field::Scalar>;

    f64 _c_mu = 0.09;
    f64 _c_eps1 = 1.44;
    f64 _c_eps2 = 1.92;
    f64 _sigma_k = 1.0;
    f64 _sigma_eps = 1.3;

    SharedPtr<field::Scalar> _k, _epsilon, _mu, _mu_t;
    SharedPtr<field::Velocity> _U, _mdot;
    SharedPtr<field::Density> _rho;
    f64 _dt = 0.0;
    bool _is_transient = false;
};

} // namespace prism::turbulence
