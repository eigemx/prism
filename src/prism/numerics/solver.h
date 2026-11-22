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

struct SolverOptions {
    // Maximum number of solver iterations (driver iterations).
    size_t max_iterations {10};

    // Convergence tolerance for the Moukalled-scaled residual
    // |Ax - b| / max(|diag(A) * x|).
    f64 tolerance {1e-7};

    // Optional reference cell used to pin the constant mode of an
    // all-Neumann (closed-domain) system.
    Optional<size_t> reference_cell {NullOption};

    // Value the reference cell is anchored to. Only used when reference_cell is set.
    f64 reference_value {0.0};

    // If true, zero the equation matrix and RHS after solving. Set false when the
    // caller still needs the assembled/relaxed matrix, e.g. to build D = V/aP or
    // to perform an explicit H-by-A update.
    bool clear_coeffs {true};

    // If true, skip preconditioner preparation and reuse the factorization from the
    // previous solve. Valid only when the matrix values are unchanged since the last
    // prepare(), e.g. successive non-orthogonal correctors of one pressure matrix.
    bool reuse_prepared {false};
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
     * @param options Iteration, tolerance, reference-cell and matrix-lifecycle options. */
    virtual auto solve(eqn::Transport& eq, const SolverOptions& options = {}) -> SolverResult;

    /** @brief Re-attaches the solver to A without recomputing the factorization.
     *
     * Used with SolverOptions::reuse_prepared: the equation matrix was re-assembled (and may
     * have been reallocated), so solvers that alias its storage must re-point to it. The
     * factorization is still valid because the matrix values are unchanged. */
    virtual void reattach(const SparseMatrix& A) { (void)A; }

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

/** @brief Solvers used by the pressure-velocity coupling algorithms.
 *
 * The algorithm owns these solvers for its lifetime, so a solver instance persists across
 * iterations and time steps. The caller is responsible for choosing solvers that match the
 * equation: momentum is nonsymmetric, pressure is symmetric positive (semi)definite.
 */
struct LinearSolvers {
    // Momentum-equation solver.
    UniquePtr<ISolver> momentum;

    // Pressure-equation solver.
    UniquePtr<ISolver> pressure;

    auto valid() const noexcept -> bool { return momentum != nullptr && pressure != nullptr; }
};

} // namespace prism::solver
