#pragma once
#include <Eigen/IterativeLinearSolvers>
#include <fstream>

#include "prism/equation/transport.h"
#include "prism/log.h"

inline void writeToCSVfile(const std::string& name, const auto& matrix) {
    std::ofstream file(name.c_str());
    file << matrix;
}

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
        -> IterationData = 0;
};

template <typename Field, typename Relaxer>
class BiCGSTAB : public ISolver<Field> {
  public:
    auto solve(eqn::Transport<Field>& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-4,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) -> IterationData override;
};

/*
template <typename Field, typename Relaxer>
class GaussSeidel : public ISolver<Field> {
  public:
    void solve(eqn::Transport<Field>& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-4,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) override;
};
*/

template <typename Field, typename Relaxer>
auto BiCGSTAB<Field, Relaxer>::solve(eqn::Transport<Field>& eqn,
                                     std::size_t n_iter,
                                     double eps,
                                     double lambda) -> IterationData {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();
    auto& phi = eqn.field();
    auto& phi_prev = eqn.prevIterField();

    auto init_res = 0.0;
    double current_res = 0.0;
    bool converged = false;
    std::size_t iter = 0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    Relaxer rx;

    for (; iter < n_iter; iter++) {
        eqn.updateCoeffs();
        rx.preRelax(eqn, lambda);

        if (iter == 0) {
            init_res = detail::residual(A, phi.values(), b);
        }

        // do at least one iteration
        phi_prev.values() = phi.values();
        phi.values() = bicg.compute(A).solveWithGuess(b, phi.values());
        rx.postRelax(eqn, lambda);

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

/*
template <typename Field, typename Relaxer>
void GaussSeidel<Field, Relaxer>::solve(eqn::Transport<Field>& eqn,
                                        std::size_t n_iter,
                                        double eps,
                                        double lambda) {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();

    auto& phi = eqn.field();
    auto& phi_prev = eqn.field_prev_iter();

    Relaxer rx;

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.updateCoeffs();
        rx.preRelax(eqn, lambda);

        // calculate the residuals and its norm
        auto res = (A * phi.values()) - b;
        auto res_norm = res.norm();

        // check for convergence
        if (res_norm < eps) {
            log::info("Converged after {} iterations", i);
            log::info("Residual: {}", res_norm);
            break;
        }

        phi_prev.values() = phi.values();

        // Implementation of Gauss-Seidel
        for (int j = 0; j < phi.values().size(); j++) {
            double sum = 0.0;
            for (int k = 0; k < phi.values().size(); k++) {
                if (k != j) {
                    sum += A.coeff(j, k) * phi.values()(k);
                }
            }
            phi.values()(j) = (b(j) - sum) / A.coeff(j, j);
        }

        rx.postRelax(eqn, lambda);

        log::info("Iteration: {}, Residual: {}", i, res_norm);

        // zero out the left & right hand side vector, for the next iteration
        eqn.zeroOutCoeffs();
    }
}
*/

} // namespace prism::solver
