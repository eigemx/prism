#pragma once

#include <memory>

#include "boundary.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/linear.h"
#include "prism/log.h"
#include "prism/scheme/convection/convection.h"
#include "prism/scheme/diffusion/diffusion.h"
#include "prism/scheme/scheme.h"
#include "prism/scheme/source/source.h"
#include "prism/scheme/temporal/temporal.h"

namespace prism::eqn {

/// TODO: Make Transport check consistincey of the conserved transport field, meaning that a
/// transport equation shall have only one transport field.

/// TODO: Schemes that require no correction shall be updated only once.

/// TODO: Transport should have at most one convection scheme, at most one diffusion scheme, and
/// at most one temporal scheme.

class Transport : public LinearSystem,
                  public prism::boundary::BHManagerProvider<
                      eqn::boundary::IEquationBoundaryHandler<Transport>> {
  public:
    template <scheme::ISchemeBased Scheme, scheme::ISchemeBased... Schemes>
    Transport(Scheme&& scheme, Schemes&&... schemes);

    void updateCoeffs();
    void zeroOutCoeffs();

    auto underRelaxFactor() const -> f64;
    void setUnderRelaxFactor(f64 factor);
    void relax();

    auto field() const -> const SharedPtr<field::Scalar>&;

    /// TODO: maybe returns a casted pointer to the scheme (convection or diffusion) instead of
    /// the full scheme?
    /// TODO: since we delegate adding schemes to overloaded addScheme for each scheme, we can
    /// directly store as SharedPtr<IAppliedDiffusion> or SharedPtr<IAppliedConvection> and avoid
    /// need for castScheme later.
    auto convectionScheme() -> SharedPtr<scheme::IFullScheme>;
    auto diffusionScheme() -> SharedPtr<scheme::IFullScheme>;
    auto temporalScheme() -> SharedPtr<scheme::IFullScheme>;

    auto isTransient() const noexcept -> bool;

    template <typename Scheme>
    void addScheme(Scheme&& scheme);

    template <scheme::convection::IAppliedConvectionBased Convection>
    void addScheme(Convection&& convection);

    template <scheme::diffusion::IAppliedDiffusionBased Diffusion>
    void addScheme(Diffusion&& diffusion);

    template <scheme::source::IExplicitSourceBased Source>
    void addScheme(Source&& source);

    template <scheme::source::IImplicitSourceBased Source>
    void addScheme(Source&& source);

    template <scheme::temporal::ITemporalBased Temporal>
    void addScheme(Temporal&& temporal);

  private:
    /// TODO: Should this be IConvectionBased and IDiffusionBased instead?
    // Conserved field of the equation
    SharedPtr<field::Scalar> _phi;

    std::vector<SharedPtr<scheme::IFullScheme>> _schemes;
    std::vector<SharedPtr<scheme::source::IExplicitSource>> _sources;

    SharedPtr<scheme::IFullScheme> _conv_scheme;
    SharedPtr<scheme::IFullScheme> _diff_scheme;
    SharedPtr<scheme::IFullScheme> _temporal_scheme;

    size_t _n_corrected_schemes {0};

    f64 _relaxation_factor {1.0};

    bool _is_transient {false};
};

class Momentum : public Transport {
  public:
    using Transport::Transport;
};

template <scheme::ISchemeBased Scheme, scheme::ISchemeBased... Schemes>
Transport::Transport(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()), LinearSystem(scheme.field()->mesh()->cellCount()) {
    /// TODO: when first scheme is temporal, eqn.field() returns zero-valued field. Check this.
    // add the first mandatory scheme
    addScheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (addScheme(std::forward<Schemes>(schemes)), ...);
}
template <typename Scheme>
void Transport::addScheme(Scheme&& scheme) {
    if (scheme.needsCorrection()) {
        _n_corrected_schemes++;
    }

    _schemes.emplace_back(std::make_shared<Scheme>(std::forward<Scheme>(scheme)));
}

template <scheme::convection::IAppliedConvectionBased Convection>
void Transport::addScheme(Convection&& convection) {
    log::debug("Transport::addScheme() found a convection scheme.");
    if (convection.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto conv_scheme = std::make_shared<Convection>(std::forward<Convection>(convection));
    _conv_scheme = conv_scheme;
    _schemes.emplace_back(conv_scheme);
}

template <scheme::diffusion::IAppliedDiffusionBased Diffusion>
void Transport::addScheme(Diffusion&& diffusion) {
    log::debug("Transport::addScheme() found a diffusion scheme.");
    if (diffusion.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto diff_scheme = std::make_shared<Diffusion>(std::forward<Diffusion>(diffusion));
    _diff_scheme = diff_scheme;
    _schemes.emplace_back(diff_scheme);
}

template <scheme::source::IExplicitSourceBased Source>
void Transport::addScheme(Source&& source) {
    log::debug("Transport::addScheme() found an explicit source scheme.");
    if (source.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto src_scheme = std::make_shared<Source>(std::forward<Source>(source));
    _sources.emplace_back(src_scheme);
}

template <scheme::source::IImplicitSourceBased Source>
void Transport::addScheme(Source&& source) {
    log::debug("Transport::addScheme() found an implicit source scheme.");
    if (source.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto src_scheme = std::make_shared<Source>(std::forward<Source>(source));
    _schemes.emplace_back(src_scheme);
}

template <scheme::temporal::ITemporalBased Temporal>
void Transport::addScheme(Temporal&& temporal) {
    log::debug("Transport::addScheme() found a temporal scheme.");
    if (temporal.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto temp_scheme = std::make_shared<Temporal>(std::forward<Temporal>(temporal));
    _schemes.emplace_back(temp_scheme);
    _is_transient = true;
}
} // namespace prism::eqn
