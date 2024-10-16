#pragma once

#include <concepts>
#include <memory>

#include "boundary.h"
#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/field/scalar.h"
#include "prism/linear.h"
#include "prism/log.h"
#include "prism/scheme/convection.h"
#include "prism/scheme/diffusion.h"
#include "prism/scheme/source.h"

namespace prism::eqn {

// TODO: Make Transport check consistincey of the  conserved transport ScalarField meaning
// that a transport equation shall have only one transport field.

// TODO: Schemes that require no correction shall be updated only once.
// TODO: Transport should have at most one convection scheme and at most one diffusion
// scheme
template <field::IScalarBased Field, typename BHManagerSetter>
class GeneralTransport {
  private:
    BHManagerSetter _setter;
};

template <field::IScalarBased Field = field::Scalar>
class Transport : public LinearSystem,
                  public prism::boundary::BHManagersProvider<
                      eqn::boundary::IEquationBoundaryHandler<Transport<Field>>> {
  public:
    using FieldType = Field;

    template <typename Scheme, typename... Schemes>
    Transport(Scheme&& scheme, Schemes&&... schemes);

    void updateCoeffs();
    void zeroOutCoeffs();

    auto inline field() const -> const Field& { return _phi; }
    auto inline field() -> Field& { return _phi; }
    auto convectionScheme() -> SharedPtr<scheme::IFullScheme<Field>> { return _conv_scheme; }

    auto diffusionScheme() -> SharedPtr<scheme::IFullScheme<Field>> { return _diff_scheme; }

    template <typename Scheme>
    void addScheme(Scheme&& scheme);

  private:
    template <typename Convection>
        requires std::derived_from<
            Convection,
            scheme::convection::IConvection<typename Convection::RhoType, Field>>
    void addScheme(Convection&& convection);

    template <typename Diffusion>
        requires std::derived_from<
            Diffusion,
            scheme::diffusion::IAppliedDiffusion<typename Diffusion::KappaType,
                                                 typename Diffusion::FieldType>>
    void addScheme(Diffusion&& diffusion);

    template <typename Source>
        requires std::derived_from<Source, scheme::source::IExplicitSource>
    void addScheme(Source&& source);


    std::vector<SharedPtr<scheme::IFullScheme<Field>>> _schemes;
    std::vector<SharedPtr<scheme::source::IExplicitSource>> _sources;
    Field _phi;     // Conserved field of the equation
    SharedPtr<scheme::IFullScheme<Field>> _conv_scheme = nullptr;
    SharedPtr<scheme::IFullScheme<Field>> _diff_scheme = nullptr;
    std::size_t _n_corrected_schemes {0};
};


template <field::IScalarBased Field>
template <typename Scheme, typename... Schemes>
Transport<Field>::Transport(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()),
      LinearSystem(scheme.field().mesh().cellCount()) {
    // add the first mandatory scheme
    addScheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (addScheme(std::forward<Schemes>(schemes)), ...);

    if (_diff_scheme) {
        log::debug("Transport::addScheme() found a diffusion scheme.");
    }

    if (_conv_scheme) {
        log::debug("Transport::addScheme() found a convection scheme.");
    }
}


template <field::IScalarBased Field>
void Transport<Field>::updateCoeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all equation's finite volume schemes
    for (auto& scheme : _schemes) {
        // apply the scheme
        scheme->apply();

        // update quation's universal coefficient matrix and RHS vector
        matrix() += scheme->matrix();
        rhs() += scheme->rhs();
    }

    // TODO: this does not consider implicit sources
    for (auto& source : _sources) {
        source->apply();
        rhs() += source->rhs();
    }

    prism::boundary::detail::applyBoundaryIfExists("prism::eqn::Transport", *this);
}


template <field::IScalarBased Field>
void Transport<Field>::zeroOutCoeffs() {
    for (auto& scheme : _schemes) {
        scheme->matrix().setZero();
        scheme->rhs().setZero();
    }

    // TODO: this does not consider implicit sources
    for (auto& source : _sources) {
        source->rhs().setZero();
    }

    // zero out the universal coefficient matrix and RHS vector
    matrix().setZero();
    rhs().setZero();
}

template <field::IScalarBased Field>
template <typename Scheme>
void Transport<Field>::addScheme(Scheme&& scheme) {
    if (scheme.needsCorrection()) {
        _n_corrected_schemes++;
    }

    _schemes.emplace_back(std::make_shared<Scheme>(std::forward<Scheme>(scheme)));
}

template <field::IScalarBased Field>
template <typename Convection>
    requires std::derived_from<
        Convection,
        scheme::convection::IConvection<typename Convection::RhoType, Field>>
void Transport<Field>::addScheme(Convection&& convection) {
    if (convection.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto conv_scheme = std::make_shared<Convection>(std::forward<Convection>(convection));
    _conv_scheme = conv_scheme;
    _schemes.emplace_back(conv_scheme);
}

template <field::IScalarBased Field>
template <typename Diffusion>
    requires std::derived_from<
        Diffusion,
        scheme::diffusion::IAppliedDiffusion<typename Diffusion::KappaType,
                                             typename Diffusion::FieldType>>
void Transport<Field>::addScheme(Diffusion&& diffusion) {
    if (diffusion.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto diff_scheme = std::make_shared<Diffusion>(std::forward<Diffusion>(diffusion));
    _diff_scheme = diff_scheme;
    _schemes.emplace_back(diff_scheme);
}

template <field::IScalarBased Field>
template <typename Source>
    requires std::derived_from<Source, scheme::source::IExplicitSource>
void Transport<Field>::addScheme(Source&& source) {
    if (source.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto src_scheme = std::make_shared<Source>(std::forward<Source>(source));
    _sources.emplace_back(src_scheme);
}

} // namespace prism::eqn
