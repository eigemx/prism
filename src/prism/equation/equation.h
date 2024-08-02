#pragma once

#include <memory>
#include <string>

#include "prism/field/field.h"
#include "prism/linear.h"
#include "prism/mesh/pmesh.h"
#include "prism/schemes/fvscheme.h"
#include "prism/types.h"

namespace prism {

// TODO: Make TransportEquation check consistincey of the  conserved transport ScalarField meaning
// that a transport equation shall have only one transport field.

// TODO: Schemes that require no correction shall be updated only once.

// TODO: Clean the below clutter!

template <typename Field = field::Scalar>
class TransportEquation : public LinearSystem {
  public:
    template <typename Scheme, typename... Schemes>
    TransportEquation(Scheme&& scheme, Schemes&&... schemes)
        : _phi(scheme.field().value()),
          _phi_old(scheme.field().value()),
          LinearSystem(scheme.field().value().mesh().n_cells()) {
        // add the first mandatory scheme
        add_scheme(std::forward<Scheme>(scheme));

        // add the rest of the schemes, if any
        (add_scheme(std::forward<Schemes>(schemes)), ...);
    }

    void update_coeffs() {
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
    void zero_out_coeffs() {
        for (auto& scheme : _schemes) {
            scheme->matrix().setZero();
            scheme->rhs().setZero();
        }

        // zero out the universal coefficient matrix and RHS vector
        matrix().setZero();
        rhs().setZero();
    }

    auto inline field() const -> const field::Scalar& { return _phi; }
    auto inline field() -> field::Scalar& { return _phi; }

    auto inline field_prev_iter() const -> const field::Scalar& { return _phi_old; }
    auto inline field_prev_iter() -> field::Scalar& { return _phi_old; }

    template <typename S>
    void add_scheme(S&& scheme) {
        if (scheme.requires_correction()) {
            _n_corrected_schemes++;
        }

        _schemes.emplace_back(std::make_shared<S>(std::forward<S>(scheme)));
    }

    template <typename S>
    void add_scheme(S& scheme) {
        if (scheme.requires_correction()) {
            _n_corrected_schemes++;
        }

        _schemes.emplace_back(std::make_shared<S>(scheme));
    }

  private:
    std::vector<std::shared_ptr<scheme::FVScheme<Field>>> _schemes;
    Field _phi;     // Conserved field of the equation
    Field _phi_old; // Previous iteration value of the field

    std::size_t _n_corrected_schemes {0};
};

} // namespace prism