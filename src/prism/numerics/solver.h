#pragma once

#include <Eigen/IterativeLinearSolvers>

#include "prism/equation/transport.h"
#include "prism/log.h"
#include "prism/numerics/relax.h"

namespace prism::solver::detail {

auto inline residual(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> double {
    auto ac_phic = A.diagonal().cwiseProduct(x);
    auto res_scaled = ((A * x) - b).cwiseAbs() / (ac_phic.maxCoeff() + EPSILON);

    return res_scaled.maxCoeff();
}

} // namespace prism::solver::detail

namespace prism::solver {

class IterationData {
  public:
    IterationData(std::size_t iteration, double initial_residual, double final_residual);
    void setAsConverged() noexcept;
    auto atConvergence() const noexcept -> bool;
    auto iteration() const noexcept -> std::size_t;
    auto initialResidual() const noexcept -> double;
    auto finalResidual() const noexcept -> double;

  private:
    std::size_t _iteration {};
    double _initial_residual {};
    double _final_residual {};
    bool _converged {false};
};

template <typename Field>
class ISolver {
  public:
    virtual auto solve(eqn::Transport<Field>& eq, std::size_t n_iter, double eps, double lambda)
        -> IterationData;
    virtual auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b)
        -> VectorXd = 0;

  private:
    ImplicitUnderRelaxation<Field> _relaxer;
};

template <typename Field>
class BiCGSTAB : public ISolver<Field> {
  public:
    auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd override;
};

template <typename Field>
class GaussSeidel : public ISolver<Field> {
  public:
    auto step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b) -> VectorXd override;
};

template <typename Field>
auto ISolver<Field>::solve(eqn::Transport<Field>& eqn,
                           std::size_t n_iter,
                           double eps,
                           double lambda) -> IterationData {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();
    auto& phi = eqn.field();

    auto init_res = 0.0;
    double current_res = 0.0;
    bool converged = false;
    std::size_t iter = 0;

    for (; iter < n_iter; iter++) {
        eqn.updateCoeffs();
        _relaxer.preRelax(eqn, lambda);

        if (iter == 0) {
            init_res = detail::residual(A, phi.values(), b);
        }

        // do at least one iteration
        phi.values() = step(A, phi.values(), b);
        current_res = detail::residual(A, phi.values(), b);

        // check for convergence
        if (current_res < eps) {
            converged = true;
            // zero out coefficients by default, in case user calls updateCoeffs() later
            eqn.zeroOutCoeffs();
            iter++;
            break;
        }

        // zero out the left & right hand side vector, for the next iteration
        eqn.zeroOutCoeffs();
    }

    log::info("Residuals: Initial = {:.4e} | Final: {:.4e} (nIterations = {})",
              init_res,
              current_res,
              iter);

    return {n_iter, init_res, current_res};
}

template <typename Field>
auto BiCGSTAB<Field>::step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b)
    -> VectorXd {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    return bicg.compute(A).solveWithGuess(b, x);
}

template <typename Field>
auto GaussSeidel<Field>::step(const SparseMatrix& A, const VectorXd& x, const VectorXd& b)
    -> VectorXd {
    VectorXd x_new = x;
    for (int j = 0; j < x.size(); j++) {
        double sum = 0.0;
        for (int k = 0; k < x.size(); k++) {
            if (k != j) {
                sum += A.coeff(j, k) * x[k];
            }
        }
        x_new[j] = (b[j] - sum) / A.coeff(j, j);
    }
    return x_new;
}

} // namespace prism::solver
