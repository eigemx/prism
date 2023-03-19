#pragma once

#include <memory>

#include "../equation.h"
#include "../field.h"
#include "../types.h"

namespace prism::solver {

class SolverBase {
  public:
    SolverBase() = default;
    SolverBase(const SolverBase& s) = delete;
    SolverBase(SolverBase&& s) = delete;
    auto operator=(SolverBase&& s) -> SolverBase& = delete;
    auto operator=(const SolverBase& s) -> SolverBase& = delete;
    virtual ~SolverBase() = default;

    virtual void solve(Equation& eq, std::size_t n_iter, double eps) = 0;
};

class GaussSeidel : public SolverBase {
  public:
    void solve(Equation& eq, std::size_t n_iter = 1000, double eps = 1e-6) override;
};

} // namespace prism::solver