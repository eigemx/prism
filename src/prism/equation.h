#pragma once

#include <string>

#include "field.h"
#include "fvscheme.h"
#include "linear.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation : public LinearSystem {
  public:
    // TODO: is it corect to call Scheme&& an rvalue reference?
    // Equation constructor acceptes only an rvalue reference to a scheme. Because each
    // equation should own its schemes.
    template <typename Scheme, typename... Schemes>
    Equation(Scheme&& scheme, Schemes&&... schemes);

    void update_coeffs();
    void zero_out_coeffs();

    auto inline field() const -> const ScalarField& { return _phi; }
    auto inline field() -> ScalarField& { return _phi; }

    // previous iteration value of the scalar field
    auto inline field_prev_iter() const -> const ScalarField& { return _phi_old; }
    auto inline field_prev_iter() -> ScalarField& { return _phi_old; }

    template <typename S>
    void add_scheme(S&& scheme) {
        const auto& field = scheme.field();

        if (&field.mesh() != &_phi.mesh()) {
            throw std::runtime_error(
                "Equation::add_scheme() was called with a scheme defined over a different mesh.");
        }

        _schemes.emplace_back(std::make_shared<S>(std::forward<S>(scheme)));
    }

    template <typename S>
    void add_scheme(S& scheme) {
        const auto& field = scheme.field();

        if (&field.mesh() != &_phi.mesh()) {
            throw std::runtime_error(
                "Equation::add_scheme() was called with a scheme defined over a different mesh.");
        }

        _schemes.emplace_back(std::make_shared<S>(scheme));
    }

  private:
    std::vector<std::shared_ptr<FVScheme>> _schemes;
    ScalarField& _phi;
    ScalarField _phi_old;
};

template <typename Scheme, typename... Schemes>
Equation::Equation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()),
      _phi_old(scheme.field()),
      LinearSystem(scheme.field().mesh().n_cells()) {
    _schemes.reserve(sizeof...(Schemes) + 1);
    add_scheme(std::forward<Scheme>(scheme));
    (add_scheme(std::forward<Schemes>(schemes)), ...);

    auto n_cells = _phi.mesh().n_cells();

    assert(n_cells > 0 &&
           "Equation constructor was called for a scalar field defined over an empty mesh.");

    matrix() = SparseMatrix(n_cells, n_cells);
    rhs() = VectorXd::Zero(n_cells);
}

} // namespace prism