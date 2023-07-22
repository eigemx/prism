#pragma once

#include <string>

#include "field.h"
#include "fvscheme.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation {
  public:
    // Equation constructor acceptes only an rvalue reference to a scheme. Because each
    // equation should own its schemes.
    template <typename Scheme, typename... Schemes>
    Equation(Scheme&& scheme, Schemes&&... schemes);

    void update_coeffs();
    void zero_out_coeffs();

    auto inline coeff_matrix() const -> const SparseMatrix& { return _coeff_matrix; }
    auto inline coeff_matrix() -> SparseMatrix& { return _coeff_matrix; }

    auto inline rhs_vector() const -> const VectorXd& { return _rhs_vector; }
    auto inline rhs_vector() -> VectorXd& { return _rhs_vector; }

    auto inline scalar_field() const -> const ScalarField& { return _phi; }
    auto inline scalar_field() -> ScalarField& { return _phi; }

    auto inline scalar_field_old() const -> const ScalarField& { return _phi_old; }
    auto inline scalar_field_old() -> ScalarField& { return _phi_old; }

    void relax(double omega);


  private:
    SparseMatrix _coeff_matrix;
    VectorXd _rhs_vector;

    std::vector<std::shared_ptr<FVScheme>> _schemes;
    ScalarField& _phi;
    ScalarField _phi_old;
};

template <typename Scheme, typename... Schemes>
Equation::Equation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()), _phi_old(scheme.field()) {
    _schemes.reserve(sizeof...(Schemes) + 1);
    _schemes.emplace_back(std::make_shared<Scheme>(std::forward<Scheme>(scheme)));
    (_schemes.emplace_back(std::make_shared<Schemes>(std::forward<Schemes>(schemes))), ...);

    auto n_cells = _phi.mesh().n_cells();

    assert(n_cells > 0 &&
           "Equation constructor was called for a scalar field defined over an empty mesh.");

    _coeff_matrix = SparseMatrix(n_cells, n_cells);
    _rhs_vector = VectorXd::Zero(n_cells);
}

} // namespace prism