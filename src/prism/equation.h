#pragma once

#include <string>

#include "field.h"
#include "fvscheme.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation {
  public:
    Equation(ScalarField& phi, std::vector<FVScheme*> schemes);

    void update_coeffs();
    void zero_out_coeffs();

    auto inline coeff_matrix() const -> const SparseMatrix& { return _coeff_matrix; }
    auto inline coeff_matrix() -> SparseMatrix& { return _coeff_matrix; }

    auto inline rhs_vector() const -> const VectorXd& { return _rhs_vector; }
    auto inline rhs_vector() -> VectorXd& { return _rhs_vector; }

    auto inline scalar_field() const -> const ScalarField& { return _phi; }
    auto inline scalar_field() -> ScalarField& { return _phi; }


  private:
    SparseMatrix _coeff_matrix;
    VectorXd _rhs_vector;

    std::vector<FVScheme*> _schemes;
    ScalarField& _phi;
};
} // namespace prism