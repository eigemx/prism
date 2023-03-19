#pragma once

#include <memory>
#include <string>

#include "diffusion/diffusion.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation {
  public:
    Equation(ScalarField& phi) : _phi(phi) {};
    Equation(ScalarField& phi, std::vector<FVScheme*> schemes);

    void update_coeffs();

    auto coeff_matrix() const -> const SparseMatrix& { return _unified_coeff_matrix; }
    auto coeff_matrix() -> SparseMatrix& { return _unified_coeff_matrix; }

    auto rhs_vector() const -> const VectorXd& { return _unified_rhs_vector; }
    auto rhs_vector() -> VectorXd& { return _unified_rhs_vector; }

    auto scalar_field() const -> const ScalarField& { return _phi; }
    auto scalar_field() -> ScalarField& { return _phi; }

    void reset();

  private:
    SparseMatrix _unified_coeff_matrix;
    VectorXd _unified_rhs_vector;
    std::vector<FVScheme*> _schemes;
    ScalarField& _phi;
};
} // namespace prism