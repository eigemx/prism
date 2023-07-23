#pragma once

#include "../equation.h"

namespace prism::solver {

class RelaxerBase {
  public:
    RelaxerBase() = default;
    RelaxerBase(double lambda) : _lambda(lambda) {}

    virtual void pre_relax(Equation& eqn) const = 0;
    virtual void post_relax(Equation& eqn) const = 0;

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
    void inline post_relax(Equation& eqn) const override {}
};


class ExplicitUnderRelaxation : public RelaxerBase {
  public:
    ExplicitUnderRelaxation() = default;
    ExplicitUnderRelaxation(double lambda = 0.9) : RelaxerBase(lambda) {}

    void inline pre_relax(Equation& eqn) const override {}
    void post_relax(Equation& eqn) const override;
};

} // namespace prism::solver