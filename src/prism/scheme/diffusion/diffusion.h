#pragma once

#include <cassert>

#include "diffusion_boundary.h"
#include "nonortho.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/utilities.h"
#include "prism/scheme/boundary.h"
#include "prism/scheme/scheme.h"
#include "prism/types.h"


namespace prism::scheme::diffusion {
// Basic base class for all diffusion schemes, without templating clutter.
class IDiffusion {};

// Base class for all diffusion schemes with shared methods (apply(), field(), and kappa()).
template <typename Kappa>
class IAppliedDiffusion : public IFullScheme {
  public:
    IAppliedDiffusion(const SharedPtr<Kappa>& kappa, const SharedPtr<field::Scalar>& phi);
    auto kappa() -> const SharedPtr<Kappa>& { return _kappa; }

    using KappaType = Kappa;

  private:
    SharedPtr<Kappa> _kappa;
};

// Concept for diffusion schemes that are based on IAppliedDiffusion.
template <typename T>
concept IAppliedDiffusionBased = std::derived_from<T, IAppliedDiffusion<typename T::KappaType>>;

// Base class for all diffusion schemes that do not need non-orthogonal correction.
class INonCorrected : public IDiffusion {};

// Base class for all diffusion schemes that need non-orthogonal correction.
class ICorrected : public IDiffusion {};

// Diffusion (laplacian) scheme without non-orthogonal correction.
template <typename Kappa = field::Scalar>
class NonCorrected : public INonCorrected,
                     public IAppliedDiffusion<Kappa>,
                     public prism::boundary::BHManagerProvider<
                         boundary::ISchemeBoundaryHandler<NonCorrected<Kappa>>> {
  public:
    NonCorrected(const SharedPtr<Kappa>& kappa, const SharedPtr<field::Scalar>& phi);
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

// Diffusion (laplacian) scheme with non-orthogonal correction.
template <typename Kappa = field::Scalar,
          typename NonOrthoCorrector = nonortho::OverRelaxedCorrector>
class Corrected : public ICorrected,
                  public IAppliedDiffusion<Kappa>,
                  public prism::boundary::BHManagerProvider<
                      boundary::ISchemeBoundaryHandler<Corrected<Kappa, NonOrthoCorrector>>> {
  public:
    Corrected(const SharedPtr<Kappa>& kappa, const SharedPtr<field::Scalar>& phi);

    auto corrector() const noexcept -> const NonOrthoCorrector& { return _corrector; }
    auto needsCorrection() const noexcept -> bool override { return true; }

    using NonOrthoCorrectorType = NonOrthoCorrector;

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;

    NonOrthoCorrector _corrector;
};

//
// diffusion::IAppliedDiffusion implementation
//
template <typename Kappa>
IAppliedDiffusion<Kappa>::IAppliedDiffusion(const SharedPtr<Kappa>& kappa,
                                            const SharedPtr<field::Scalar>& phi)
    : _kappa(kappa), IFullScheme(phi) {}

///
/// diffusion::NonCorrected implementation
///
template <typename Kappa>
NonCorrected<Kappa>::NonCorrected(const SharedPtr<Kappa>& kappa,
                                  const SharedPtr<field::Scalar>& phi)
    : IAppliedDiffusion<Kappa>(kappa, phi) {
    // add default boundary handlers for NonCorrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <typename Kappa>
void NonCorrected<Kappa>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::NonCorrected", *this);
}

template <typename Kappa>
void NonCorrected<Kappa>::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field()->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.areaVector();

    // The following is based on equation (8.93) from Moukalled et al. (2015).
    // this handles the general case when kappa is a field::Tensor.
    Vector3d Sf_prime = detail::valueAtFace(this->kappa(), face) * Sf;

    // geometric diffusion coefficient
    // Taking the norm of Sf_prime discards the sign of kappa, so we use the following instead
    const f64 g_diff = Sf_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    // since we are applying the discretized equation per face basis, not per cell basis,
    // we need to insert the coefficients for both owner and neighbor cells when roles are
    // reversed.
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);
}

//
// diffusion::Corrected implementation
//
template <typename Kappa, typename NonOrthoCorrector>
Corrected<Kappa, NonOrthoCorrector>::Corrected(const SharedPtr<Kappa>& kappa,
                                               const SharedPtr<field::Scalar>& phi)
    : IAppliedDiffusion<Kappa>(kappa, phi) {
    // add default boundary handlers for Corrected
    using Scheme = Corrected<Kappa, NonOrthoCorrector>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <typename Kappa, typename NonOrthoCorrector>
void Corrected<Kappa, NonOrthoCorrector>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::Corrected", *this);
}

template <typename Kappa, typename NonOrthoCorrector>
void Corrected<Kappa, NonOrthoCorrector>::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field()->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const f64 d_CF_norm = d_CF.norm();
    const Vector3d e = d_CF / d_CF_norm;

    const auto Sf = mesh::outwardAreaVector(face, owner);
    const Vector3d Sf_prime = detail::valueAtFace(this->kappa(), face) * Sf;
    const auto [Ef_prime, Tf_prime] = _corrector.decompose(Sf_prime, e);

    // geometric diffusion coefficient
    const f64 g_diff = Ef_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

    /// g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);

    // update right hand side
    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukalled et al., 2015)
    const Vector3d grad_f = this->field()->gradAtFace(face);
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}


} // namespace prism::scheme::diffusion
