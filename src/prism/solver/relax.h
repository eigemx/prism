#pragma once

#include "../equation.h"

namespace prism::solver {

class ImplicitUnderRelaxation {
  public:
    ImplicitUnderRelaxation(double lambda = 0.9) : _lambda(lambda) {}
    void relax(Equation& eqn) const;

    auto lambda() const -> double { return _lambda; }
    auto lambda() -> double& { return _lambda; }

  private:
    double _lambda = 1.0;
};

} // namespace prism::solver