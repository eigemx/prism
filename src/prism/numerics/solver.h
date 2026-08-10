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
    auto iteration() const noexcept -> size_t;
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

    /** @brief Assembles, relaxes and solves the equation's linear system.
     * @param eq The equation to solve.
     * @param n_iter Maximum number of solver iterations.
     * @param eps Convergence tolerance (scaled residual).
     * @param ref_cell Optional reference cell used to pin the solution of singular systems.
     * @param ref_value Value to anchor the reference cell to.
     * @param clear_coeffs If true, zeroes the equation's matrix and rhs after solving. Pass
     *        false to keep the relaxed matrix, e.g. so the algorithm can reuse its diagonal to
     *        build the H-by-A coefficient tensor without re-assembling the equation. */
    virtual auto solve(eqn::Transport& eq,
                       size_t n_iter = 10,
                       f64 eps = 1e-7,
                       Optional<size_t> ref_cell = NullOption,
                       f64 ref_value = 0.0,
                       bool clear_coeffs = true) -> SolverResult;

    /** @brief Factorizes/prepares the preconditioner for the matrix A.
     *
     * Must be called once before iterating. The matrix must remain alive for the
     * duration of the iteration loop (the solver keeps a reference to it).
     *
     * @param A The matrix to factorize. */
    virtual void prepare(const SparseMatrix& A) = 0;

    /** @brief Performs one solver application using the factorization from prepare().
     *
     * For BiCGSTAB this runs one full internal Krylov solve with the cached
     * preconditioner, so the factorization performed by prepare() is reused instead of
     * being recomputed on every step.
     *
     * @param x Current solution guess.
     * @param b Right-hand side vector.
     * @return The updated solution. */
    virtual auto step(const VectorXd& x, const VectorXd& b) -> VectorXd = 0;
};

} // namespace prism::solver
