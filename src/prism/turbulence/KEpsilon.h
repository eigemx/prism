#pragma once

#include "model.h"
#include "prism/equation/transport.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::turbulence {

class IncompressibleKEpsilon : public IRASModel {
  public:
    IncompressibleKEpsilon(const SharedPtr<field::Velocity>& mdot,
                           const SharedPtr<field::Velocity>& U,
                           const SharedPtr<field::Scalar>& mu,
                           Optional<f64> dt = NullOption);

    auto k() -> const SharedPtr<field::Scalar>&;
    auto epsilon() -> const SharedPtr<field::Scalar>&;
    auto turbulentViscosity() -> const SharedPtr<field::Scalar>&;

    auto equations() -> Pair<eqn::Transport, eqn::Transport>;

  private:
    auto buildKEqn(const SharedPtr<field::Scalar>& pk) -> eqn::Transport;
    auto buildEpsilonEqn(const SharedPtr<field::Scalar>& pk) -> eqn::Transport;

    void setTurbulentViscosity();

    f64 _c_mu = 0.09;
    f64 _c_eps1 = 1.44;
    f64 _c_eps2 = 1.92;
    f64 _sigma_k = 1.0;
    f64 _sigma_eps = 1.3;

    SharedPtr<field::Scalar> _k, _epsilon, _nu, _nu_t;
    SharedPtr<field::Velocity> _U, _mdot;
    f64 _dt = 0.0;
    bool _is_transient = false;
};

} // namespace prism::turbulence
