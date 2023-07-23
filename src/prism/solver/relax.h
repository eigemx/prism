#pragma once

#include "../equation.h"

namespace prism::solver {

class RelaxerBase {
  public:
    RelaxerBase() = default;
    RelaxerBase(double lambda) : _lambda(lambda) {}

    // pre_relax() is called before solving the equation (suitable for implicit schemes)
    virtual void pre_relax(Equation& eqn) const = 0;

    // post_relax() is called after solving the equation (suitable for explicit schemes)
    virtual void post_relax(Equation& eqn) const = 0;

    // getter and setter for lambda (the under-relaxation factor)
    auto lambda() const -> double { return _lambda; }
    auto lambda() -> double& { return _lambda; }

  private:
    double _lambda = 1.0;
};


// TODO: Test this, under-relaxation when applied to basic cases of diffusion or advection
// yields wrong results. (check duct and torus cases)
class ImplicitUnderRelaxation : public RelaxerBase {
  public:
    ImplicitUnderRelaxation() = default;
    ImplicitUnderRelaxation(double lambda = 0.9) : RelaxerBase(lambda) {}

    void pre_relax(Equation& eqn) const override;

    // no post_relax() is needed for implicit under-relaxation
    void inline post_relax(Equation& eqn) const override {}
};


class ExplicitUnderRelaxation : public RelaxerBase {
  public:
    ExplicitUnderRelaxation() = default;
    ExplicitUnderRelaxation(double lambda = 0.9) : RelaxerBase(lambda) {}

    // no pre_relax() is needed for explicit under-relaxation
    void inline pre_relax(Equation& eqn) const override {}

    void post_relax(Equation& eqn) const override;
};

} // namespace prism::solver