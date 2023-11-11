#pragma once

#include <memory>
#include <string>

#include "field.h"
#include "fvscheme.h"
#include "linear.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

// TODO: Make TransportEquation check of a consisten conserved transport ScalarField
// meaning that a transport equation shall have only one transport field.
// Clean the below clutter!

class TransportEquation : public LinearSystem {
  public:
    template <typename Scheme, typename... Schemes>
    TransportEquation(Scheme&& scheme, Schemes&&... schemes);

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
    ScalarField _phi;     // Conserved field of the equation
    ScalarField _phi_old; // Previous iteration value of the field

    std::size_t _n_corrected_schemes {0};
};

template <typename Scheme, typename... Schemes>
TransportEquation::TransportEquation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field().value()),
      _phi_old(scheme.field().value()),
      LinearSystem(scheme.field().value().mesh().n_cells()) {
    // add the first mandatory scheme
    add_scheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (add_scheme(std::forward<Schemes>(schemes)), ...);
}

template <typename S>
void TransportEquation::add_scheme(S&& scheme) {
    if (scheme.requires_correction()) {
        _n_corrected_schemes++;
    }

    _schemes.emplace_back(std::make_shared<S>(std::forward<S>(scheme)));
}

template <typename S>
void TransportEquation::add_scheme(S& scheme) {
    if (scheme.requires_correction()) {
        _n_corrected_schemes++;
    }

    _schemes.emplace_back(std::make_shared<S>(scheme));
}

} // namespace prism