#pragma once

#include <memory>
#include <string>

#include "diffusion/diffusion.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation {
  public:
    Equation(ScalarField& phi, std::vector<FVScheme*> schemes);

    void update_coeffs();
    void reset_coeffs();

    auto coeff_matrix() const -> const SparseMatrix& { return _coeff_matrix; }
    auto coeff_matrix() -> SparseMatrix& { return _coeff_matrix; }

    auto rhs_vector() const -> const VectorXd& { return _rhs_vector; }
    auto rhs_vector() -> VectorXd& { return _rhs_vector; }

    auto scalar_field() const -> const ScalarField& { return _phi; }
    auto scalar_field() -> ScalarField& { return _phi; }


  private:
    SparseMatrix _coeff_matrix;
    VectorXd _rhs_vector;
    std::vector<FVScheme*> _schemes;
    ScalarField& _phi;
};
} // namespace prism