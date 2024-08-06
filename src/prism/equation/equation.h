#pragma once

#include <memory>

#include "prism/field/scalar.h"
#include "prism/linear.h"
#include "prism/schemes/fvscheme.h"

namespace prism {

// TODO: Make TransportEquation check consistincey of the  conserved transport ScalarField meaning
// that a transport equation shall have only one transport field.

// TODO: Schemes that require no correction shall be updated only once.

// TODO: Clean the below clutter!

template <typename Field = field::Scalar>
class TransportEquation : public LinearSystem {
  public:
    template <typename Scheme, typename... Schemes>
    TransportEquation(Scheme&& scheme, Schemes&&... schemes);

    void updateCoeffs();
    void zeroOutCoeffs();
    auto inline field() const -> const Field& { return _phi; }
    auto inline field() -> Field& { return _phi; }

    auto inline prevIterField() const -> const Field& { return _phi_old; }
    auto inline prevIterField() -> Field& { return _phi_old; }

    template <typename Scheme>
    void addScheme(Scheme&& scheme);

  private:
    std::vector<std::shared_ptr<scheme::FVScheme<Field>>> _schemes;
    Field _phi;     // Conserved field of the equation
    Field _phi_old; // Previous iteration value of the field

    std::size_t _n_corrected_schemes {0};
};


template <typename Field>
template <typename Scheme, typename... Schemes>
TransportEquation<Field>::TransportEquation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field().value()),
      _phi_old(scheme.field().value()),
      LinearSystem(scheme.field().value().mesh().n_cells()) {
    // add the first mandatory scheme
    addScheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (addScheme(std::forward<Schemes>(schemes)), ...);
}

template <typename Field>
void TransportEquation<Field>::updateCoeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all equation's finite volume schemes
    for (auto& scheme : _schemes) {
        // apply the scheme
        scheme->apply();

        // update quation's universal coefficient matrix and RHS vector
        matrix() += scheme->matrix();
        rhs() += scheme->rhs();
    }
}

template <typename Field>
void TransportEquation<Field>::zeroOutCoeffs() {
    for (auto& scheme : _schemes) {
        scheme->matrix().setZero();
        scheme->rhs().setZero();
    }

    // zero out the universal coefficient matrix and RHS vector
    matrix().setZero();
    rhs().setZero();
}

template <typename Field>
template <typename Scheme>
void TransportEquation<Field>::addScheme(Scheme&& scheme) {
    if (scheme.needsCorrection()) {
        _n_corrected_schemes++;
    }

    _schemes.emplace_back(std::make_shared<Scheme>(std::forward<Scheme>(scheme)));
}

} // namespace prism