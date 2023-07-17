#pragma once

#include <memory>

#include "../equation.h"
#include "../field.h"
#include "../types.h"

namespace prism::solver {

// TODO: make this a template class, with the template parameter being the solver type
// TODO: solve() should return each time an iteraion is finished (with the current residual),
// until convergence. And any consecutive call to solve() after convergence shall have no effect
class SolverBase {
  public:
    SolverBase() = default;
    SolverBase(const SolverBase& s) = delete;
    SolverBase(SolverBase&& s) = delete;
    auto operator=(SolverBase&& s) -> SolverBase& = delete;
    auto operator=(const SolverBase& s) -> SolverBase& = delete;
    virtual ~SolverBase() = default;

    virtual void solve(Equation& eq, std::size_t n_iter, double eps, double lambda) = 0;
};

class BiCGSTAB : public SolverBase {
  public:
    void solve(Equation& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-6,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) override;
};

} // namespace prism::solver