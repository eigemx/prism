#pragma once
#include <Eigen/IterativeLinearSolvers>
#include <fstream>

#include "prism/equation/transport.h"
#include "prism/log.h"
#include "prism/types.h"

inline void writeToCSVfile(const std::string& name, const auto& matrix) {
    std::ofstream file(name.c_str());
    file << matrix;
}

namespace prism::solver {

struct IterationStepData {
    std::size_t iteration {};
    double residual {};
};

// TODO: solve() should return each time an iteraion is finished (with the current residual),
// until convergence. And any consecutive call to solve() after convergence shall have no effect
template <typename Field>
class ISolver {
  public:
    virtual void solve(eqn::Transport<Field>& eq,
                       std::size_t n_iter,
                       double eps,
                       double lambda) = 0;
    // virtual auto step(Equation& eq, double eps, double lambda) -> IterationStepData = 0;
};

template <typename Field, typename Relaxer>
class BiCGSTAB : public ISolver<Field> {
  public:
    void solve(eqn::Transport<Field>& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-4,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) override;
};

template <typename Field, typename Relaxer>
class GaussSeidel : public ISolver<Field> {
  public:
    void solve(eqn::Transport<Field>& eq,
               std::size_t n_iter = 1000, // number of iterations
               double eps = 1e-4,         // convergence criteria
               double lambda = 1.0        // under-relaxation factor
               ) override;
};

template <typename Field, typename Relaxer>
void BiCGSTAB<Field, Relaxer>::solve(eqn::Transport<Field>& eqn,
                                     std::size_t n_iter,
                                     double eps,
                                     double lambda) {
    const auto& A = eqn.matrix();
    const auto& b = eqn.rhs();

    auto& phi = eqn.field();
    //auto& phi_prev = eqn.prevIterField();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    Relaxer rx;

    for (std::size_t i = 0; i < n_iter; i++) {
        eqn.updateCoeffs();
        rx.preRelax(eqn, lambda);

        // calculate the residuals and its norm
        auto res = (A * phi.values()) - b;
        auto res_norm = res.norm();

        //phi_prev.values() = phi.values();
        phi.values() = bicg.compute(A).solveWithGuess(b, phi.values());

        // check for convergence
        /*
        if (res_norm < eps) {
            log::info("Converged after {} iterations, Residual: {}", i, res_norm);
            break;
        }
        */

        rx.postRelax(eqn, lambda);

        log::info("Iteration: {}, Residual: {}", i, res_norm);

        // zero out the left & right hand side vector, for the next iteration
        eqn.zeroOutCoeffs();

    }
}

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

} // namespace prism::solver
