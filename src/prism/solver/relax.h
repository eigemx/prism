#pragma once

#include "../equation.h"

namespace prism::solver {

class RelaxerBase {
  public:
    // pre_relax() is called before solving the equation (suitable for implicit schemes)
    virtual void pre_relax(Equation& eqn, double lambda) const = 0;

    // post_relax() is called after solving the equation (suitable for explicit schemes)
    virtual void post_relax(Equation& eqn, double lambda) const = 0;
};


// TODO: Test this, under-relaxation when applied to basic cases of diffusion or advection
// yields wrong results. (check duct and torus cases)
class ImplicitUnderRelaxation : public RelaxerBase {
  public:
    void pre_relax(Equation& eqn, double lambda) const override;

    // no post_relax() is needed for implicit under-relaxation
    void inline post_relax(Equation& eqn, double lambda) const override {}
};


class ExplicitUnderRelaxation : public RelaxerBase {
  public:
    // no pre_relax() is needed for explicit under-relaxation
    void inline pre_relax(Equation& eqn, double lambda) const override {}

    void post_relax(Equation& eqn, double lambda) const override;
};

} // namespace prism::solver