#pragma once

#include <concepts>
#include <memory>

#include "prism/field/scalar.h"
#include "prism/linear.h"
#include "prism/schemes/convection.h"
#include "prism/schemes/diffusion.h"

namespace prism {

// TODO: Make TransportEquation check consistincey of the  conserved transport ScalarField meaning
// that a transport equation shall have only one transport field.

// TODO: Schemes that require no correction shall be updated only once.
// TODO: TransportEquation should have at most one convection scheme and at most one diffusion
// scheme
// TODO: TransportEquation should give read/write access to user for all its schemes

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

    auto convectionScheme() -> std::shared_ptr<scheme::FVScheme<Field>> { return _conv_scheme; }
    auto diffusionScheme() -> std::shared_ptr<scheme::FVScheme<Field>> { return _diff_scheme; }

    template <typename Scheme>
    void addScheme(Scheme&& scheme);

  private:
    template <typename Convection>
    requires std::derived_from<Convection, scheme::convection::IConvection>
    void addScheme(Convection&& convection);

    template <typename Diffusion>
    requires std::derived_from<Diffusion,
                               scheme::diffusion::IDiffusion<typename Diffusion::KappaType,
                                                             typename Diffusion::FieldType>>
    void addScheme(Diffusion&& diffusion);

    std::vector<std::shared_ptr<scheme::FVScheme<Field>>> _schemes;
    Field _phi;     // Conserved field of the equation
    Field _phi_old; // Previous iteration value of the field

    std::shared_ptr<scheme::FVScheme<Field>> _conv_scheme = nullptr;
    std::shared_ptr<scheme::FVScheme<Field>> _diff_scheme = nullptr;

    std::size_t _n_corrected_schemes {0};
};


template <typename Field>
template <typename Scheme, typename... Schemes>
TransportEquation<Field>::TransportEquation(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field().value()),
      _phi_old(scheme.field().value()),
      LinearSystem(scheme.field().value().mesh().nCells()) {
    // add the first mandatory scheme
    addScheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (addScheme(std::forward<Schemes>(schemes)), ...);

    if (_diff_scheme) {
        spdlog::debug("TransportEquation::addScheme() found a diffusion scheme.");
    }

    if (_conv_scheme) {
        spdlog::debug("TransportEquation::addScheme() found a convection scheme.");
    }
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

template <typename Field>
template <typename Convection>
requires std::derived_from<Convection, scheme::convection::IConvection>
void TransportEquation<Field>::addScheme(Convection&& convection) {
    if (convection.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto conv_scheme = std::make_shared<Convection>(std::forward<Convection>(convection));
    _conv_scheme = conv_scheme;
    _schemes.emplace_back(conv_scheme);
}

template <typename Field>
template <typename Diffusion>
requires std::derived_from<
    Diffusion,
    scheme::diffusion::IDiffusion<typename Diffusion::KappaType, typename Diffusion::FieldType>>
void TransportEquation<Field>::addScheme(Diffusion&& diffusion) {
    if (diffusion.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto diff_scheme = std::make_shared<Diffusion>(std::forward<Diffusion>(diffusion));
    _diff_scheme = diff_scheme;
    _schemes.emplace_back(diff_scheme);
}

} // namespace prism