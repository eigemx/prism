#pragma once

#include "prism/equation/transport.h"
#include "prism/types.h"

namespace prism::solver {

class SolverResult {
  public:
    SolverResult() = default;
    SolverResult(size_t iteration, f64 initial_residual, f64 final_residual);
    void setAsConverged() noexcept;
    auto hasConverged() const noexcept -> bool;
    auto iteration() const noexcept -> std::size_t;
    auto initialResidual() const noexcept -> f64;
    auto finalResidual() const noexcept -> f64;

  private:
    size_t _iteration {0};
    f64 _initial_residual {0.0};
    f64 _final_residual {0.0};
    bool _converged {false};
};

class ISolver {
  public:
    ISolver() = default;
    ISolver(const ISolver&) = default;
    ISolver(ISolver&&) = default;
    auto operator=(const ISolver&) -> ISolver& = default;
    auto operator=(ISolver&&) -> ISolver& = default;
    virtual ~ISolver() = default;

    virtual auto solve(eqn::Transport& eq, size_t n_iter = 10, f64 eps = 1e-7) -> SolverResult;
    virtual auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd = 0;
};

} // namespace prism::solver