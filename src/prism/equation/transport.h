#pragma once

#include <memory>

#include "boundary.h"
#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/field/scalar.h"
#include "prism/linear.h"
#include "prism/log.h"
#include "prism/scheme/convection/convection.h"
#include "prism/scheme/diffusion/diffusion.h"
#include "prism/scheme/source/source.h"

namespace prism::eqn {

/// TODO: Make Transport check consistincey of the conserved transport field, meaning that a
/// transport equation shall have only one transport field.

/// TODO: Schemes that require no correction shall be updated only once.

/// TODO: Transport should have at most one convection scheme and at most one diffusion scheme

template <field::IScalarBased Field = field::Scalar>
class Transport : public LinearSystem,
                  public prism::boundary::BHManagerProvider<
                      eqn::boundary::IEquationBoundaryHandler<Transport<Field>>> {
  public:
    using FieldType = Field;

    template <typename Scheme, typename... Schemes>
    Transport(Scheme&& scheme, Schemes&&... schemes);

    void updateCoeffs();
    void zeroOutCoeffs();

    auto underRelaxFactor() const -> double;
    void setUnderRelaxFactor(double factor);
    void relax();

    auto field() const -> const Field&;
    auto field() -> Field&;

    /// TODO: maybe returns a casted pointer to the scheme (convection or diffusion) instead of
    /// the full scheme?
    auto convectionScheme() -> SharedPtr<scheme::IFullScheme<Field>>;
    auto diffusionScheme() -> SharedPtr<scheme::IFullScheme<Field>>;

    template <typename Scheme>
    void addScheme(Scheme&& scheme);

    /// TODO: implement the following, for better SIMPLE and SIMPLE-like algorithms implementation
    auto AbyV() const -> field::Scalar;
    auto H() const -> field::Scalar;
    auto H1() const -> field::Scalar;

  private:
    /// TODO: Should this be IConvectionBased and IDiffusionBased instead?
    template <scheme::convection::IAppliedConvectionBased Convection>
    void addScheme(Convection&& convection);

    template <scheme::diffusion::IAppliedDiffusionBased Diffusion>
    void addScheme(Diffusion&& diffusion);

    template <scheme::source::IExplicitSourceBased Source>
    void addScheme(Source&& source);

    template <scheme::source::IImplicitSourceBased Source>
    void addScheme(Source&& source);

    std::vector<SharedPtr<scheme::IFullScheme<Field>>> _schemes;
    std::vector<SharedPtr<scheme::source::IExplicitSource>> _sources;
    Field _phi; // Conserved field of the equation
    SharedPtr<scheme::IFullScheme<Field>> _conv_scheme = nullptr;
    SharedPtr<scheme::IFullScheme<Field>> _diff_scheme = nullptr;
    std::size_t _n_corrected_schemes {0};
    double _relaxation_factor {1.0};
};


template <field::IScalarBased Field>
template <typename Scheme, typename... Schemes>
Transport<Field>::Transport(Scheme&& scheme, Schemes&&... schemes)
    : _phi(scheme.field()), LinearSystem(scheme.field().mesh()->cellCount()) {
    // add the first mandatory scheme
    addScheme(std::forward<Scheme>(scheme));

    // add the rest of the schemes, if any
    (addScheme(std::forward<Schemes>(schemes)), ...);
}


template <field::IScalarBased Field>
void Transport<Field>::updateCoeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all equation's finite volume schemes
    for (auto& scheme : _schemes) {
        // apply the scheme
        scheme->apply();

        // update equation's universal coefficient matrix and RHS vector
        matrix() += scheme->matrix();
        rhs() += scheme->rhs();
    }

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

    /// TODO: this does not consider implicit sources
    for (auto& source : _sources) {
        source->rhs().setZero();
    }

    // zero out the universal coefficient matrix and RHS vector
    matrix().setZero();
    rhs().setZero();
}

template <field::IScalarBased Field>
auto Transport<Field>::underRelaxFactor() const -> double {
    return _relaxation_factor;
}

template <field::IScalarBased Field>
void Transport<Field>::setUnderRelaxFactor(double factor) {
    if (factor < 0.0 || factor > 1.0) {
        log::warn(
            "Transport::setUnderRelaxFactor(): relaxation factor must be in [0, 1] range. "
            "Ignoring the value {} and using 1.0 instead.",
            factor);
        factor = 1.0;
        return;
    }
    _relaxation_factor = factor;
}

template <field::IScalarBased Field>
void Transport<Field>::relax() {
    auto& A = matrix();
    auto& b = rhs();
    const auto& phi = field().values();

    // Moukalled et. al, 14.2 Under-Relaxation of the Algebraic Equations
    // equation 14.9
    log::debug("Transport::relax(): applying implicit under-relaxation with factor = {}",
               _relaxation_factor);
    b += ((1.0 - _relaxation_factor) / _relaxation_factor) * A.diagonal().cwiseProduct(phi);
    A.diagonal() /= _relaxation_factor;
}

template <field::IScalarBased Field>
auto Transport<Field>::field() const -> const Field& {
    return _phi;
}

template <field::IScalarBased Field>
auto Transport<Field>::field() -> Field& {
    return _phi;
}

template <field::IScalarBased Field>
auto Transport<Field>::convectionScheme() -> SharedPtr<scheme::IFullScheme<Field>> {
    return _conv_scheme;
}

template <field::IScalarBased Field>
auto Transport<Field>::diffusionScheme() -> SharedPtr<scheme::IFullScheme<Field>> {
    return _diff_scheme;
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
template <scheme::convection::IAppliedConvectionBased Convection>
void Transport<Field>::addScheme(Convection&& convection) {
    log::debug("Transport::addScheme() found a convection scheme.");
    if (convection.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto conv_scheme = std::make_shared<Convection>(std::forward<Convection>(convection));
    _conv_scheme = conv_scheme;
    _schemes.emplace_back(conv_scheme);
}

template <field::IScalarBased Field>
template <scheme::diffusion::IAppliedDiffusionBased Diffusion>
void Transport<Field>::addScheme(Diffusion&& diffusion) {
    log::debug("Transport::addScheme() found a diffusion scheme.");
    if (diffusion.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto diff_scheme = std::make_shared<Diffusion>(std::forward<Diffusion>(diffusion));
    _diff_scheme = diff_scheme;
    _schemes.emplace_back(diff_scheme);
}

template <field::IScalarBased Field>
template <scheme::source::IExplicitSourceBased Source>
void Transport<Field>::addScheme(Source&& source) {
    log::debug("Transport::addScheme() found an explicit source scheme.");
    if (source.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto src_scheme = std::make_shared<Source>(std::forward<Source>(source));
    _sources.emplace_back(src_scheme);
}

template <field::IScalarBased Field>
template <scheme::source::IImplicitSourceBased Source>
void Transport<Field>::addScheme(Source&& source) {
    log::debug("Transport::addScheme() found an implicit source scheme.");
    if (source.needsCorrection()) {
        _n_corrected_schemes++;
    }

    auto src_scheme = std::make_shared<Source>(std::forward<Source>(source));
    _schemes.emplace_back(src_scheme);
}

} // namespace prism::eqn
