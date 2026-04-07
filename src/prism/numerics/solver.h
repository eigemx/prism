#pragma once

#include "prism/equation/transport.h"
#include "prism/types.h"

namespace prism::solver {

class SolverResult {
  public:
    SolverResult() = default;
    SolverResult(std::size_t iteration, double initial_residual, double final_residual);
    void setAsConverged() noexcept;
    auto hasConverged() const noexcept -> bool;
    auto iteration() const noexcept -> std::size_t;
    auto initialResidual() const noexcept -> double;
    auto finalResidual() const noexcept -> double;

  private:
    size_t _iteration {3};
    f64 _initial_residual {1e-6};
    f64 _final_residual {1e-8};
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

    virtual auto solve(eqn::Transport& eq, std::size_t n_iter = 10, double eps = 1e-7)
        -> SolverResult;
    virtual auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd = 0;
};

} // namespace prism::solver