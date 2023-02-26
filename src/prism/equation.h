#pragma once

#include <memory>
#include <string>

#include "diffusion/diffusion.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class LinearSystem {
  public:
    virtual auto coeff_matrix() const -> const SparseMatrix& = 0;
    virtual auto coeff_matrix() -> SparseMatrix& = 0;

    virtual auto lhs_vector() const -> const VectorXd& = 0;
    virtual auto lhs_vector() -> VectorXd& = 0;
};

class ConservedScalarSteadyEquation : public LinearSystem {
  public:
    ConservedScalarSteadyEquation(std::string scalar_name, mesh::PMesh& mesh);

    auto coeff_matrix() const -> const SparseMatrix& override { return _coeff_matrix; }
    auto coeff_matrix() -> SparseMatrix& override { return _coeff_matrix; }

    auto lhs_vector() const -> const VectorXd& override { return _b; }
    auto lhs_vector() -> VectorXd& override { return _b; }

  private:
    mesh::PMesh& _mesh;
    std::string _scalar_name;

    SparseMatrix _coeff_matrix;
    VectorXd _b;
};
} // namespace prism