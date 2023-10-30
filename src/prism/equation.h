#pragma once

#include <memory>
#include <string>

#include "field.h"
#include "fvscheme.h"
#include "linear.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class Equation : public LinearSystem {
  public:
    template <typename Scheme, typename... Schemes>
    Equation(Scheme&& scheme, Schemes&&... schemes);

    void update_coeffs();
    void zero_out_coeffs();

    auto inline field() const -> const ScalarField& { return _phi; }
    auto inline field() -> ScalarField& { return _phi; }

    auto inline field_prev_iter() const -> const ScalarField& { return _phi_old; }
    auto inline field_prev_iter() -> ScalarField& { return _phi_old; }

    template <typename S>
    void add_scheme(S&& scheme);

    template <typename S>
    void add_scheme(S& scheme);

  private:
    std::vector<std::shared_ptr<FVScheme>> _schemes;
    ScalarField& _phi;    // Conserved field of the equation
    ScalarField _phi_old; // Previous iteration value of the field

    std::size_t _n_corrected_schemes {0};
};

template <typename Scheme, typename... Schemes>
Equation::Equation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()),
      _phi_old(scheme.field()),
      LinearSystem(scheme.field().mesh().n_cells()) {
    // add the first mandatory scheme
    add_scheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (add_scheme(std::forward<Schemes>(schemes)), ...);
}

template <typename S>
void Equation::add_scheme(S&& scheme) {
    const auto& field = scheme.field();

    if (scheme.requires_correction()) {
        _n_corrected_schemes++;
    }

    if (&field.mesh() != &_phi.mesh()) {
        throw std::runtime_error(
            "Equation::add_scheme() was called with a scheme defined over a different mesh.");
    }

    _schemes.emplace_back(std::make_shared<S>(std::forward<S>(scheme)));
}

template <typename S>
void Equation::add_scheme(S& scheme) {
    const auto& field = scheme.field();

    if (scheme.requires_correction()) {
        _n_corrected_schemes++;
    }

    if (&field.mesh() != &_phi.mesh()) {
        throw std::runtime_error(
            "Equation::add_scheme() was called with a scheme defined over a different mesh.");
    }

    _schemes.emplace_back(std::make_shared<S>(scheme));
}

} // namespace prism