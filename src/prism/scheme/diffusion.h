#pragma once

#include <cassert>

#include "boundary.h"
#include "diffusion_boundary.h"
#include "nonortho.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"
#include "scheme.h"

namespace prism::scheme::diffusion {
// Basic base class for all diffusion schemes, without templating clutter.
class IDiffusion {};

// Base class for all diffusion schemes with shared methods (apply(), field(), and kappa()).
template <typename Kappa, typename Field>
class IAppliedDiffusion : public IFullScheme<Field> {
  public:
    IAppliedDiffusion(Kappa kappa, Field phi);
    auto field() -> Field override { return _phi; }
    auto kappa() -> Kappa { return _kappa; }

    using FieldType = Field;
    using KappaType = Kappa;

  private:
    Field _phi;
    Kappa _kappa;
};

// Concept for diffusion schemes that are based on IAppliedDiffusion.
template <typename T>
concept IAppliedDiffusionBased =
    std::derived_from<T, IAppliedDiffusion<typename T::KappaType, typename T::FieldType>>;

// Base class for all diffusion schemes that do not need non-orthogonal correction.
class INonCorrected : public IDiffusion {};

// Base class for all diffusion schemes that need non-orthogonal correction.
class ICorrected : public IDiffusion {};

// Diffusion (laplacian) scheme without non-orthogonal correction.
template <typename KappaType = field::UniformScalar, typename Field = field::Scalar>
class NonCorrected : public INonCorrected,
                     public IAppliedDiffusion<KappaType, Field>,
                     public prism::boundary::BHManagerProvider<
                         boundary::ISchemeBoundaryHandler<NonCorrected<KappaType, Field>>> {
  public:
    NonCorrected(KappaType kappa, Field phi);
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

// Diffusion (laplacian) scheme with non-orthogonal correction.
template <typename KappaType = field::UniformScalar,
          typename NonOrthoCorrector = nonortho::OverRelaxedCorrector,
          typename Field = field::Scalar>
class Corrected
    : public ICorrected,
      public IAppliedDiffusion<KappaType, Field>,
      public prism::boundary::BHManagerProvider<
          boundary::ISchemeBoundaryHandler<Corrected<KappaType, NonOrthoCorrector, Field>>> {
  public:
    Corrected(KappaType kappa, Field phi);

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
template <typename KappaType, typename Field>
IAppliedDiffusion<KappaType, Field>::IAppliedDiffusion(KappaType kappa, Field phi)
    : _kappa(kappa), _phi(phi), IFullScheme<Field>(phi.mesh()->cellCount()) {}

///
/// diffusion::NonCorrected implementation
///
template <typename KappaType, typename Field>
NonCorrected<KappaType, Field>::NonCorrected(KappaType kappa, Field phi)
    : IAppliedDiffusion<KappaType, Field>(kappa, phi) {
    // add default boundary handlers for NonCorrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <typename KappaType, typename Field>
void NonCorrected<KappaType, Field>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::NonCorrected", *this);
}

template <typename KappaType, typename Field>
void NonCorrected<KappaType, Field>::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field().mesh();
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
    const double g_diff = Sf_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

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
template <typename KappaType, typename NonOrthoCorrector, typename Field>
Corrected<KappaType, NonOrthoCorrector, Field>::Corrected(KappaType kappa, Field phi)
    : IAppliedDiffusion<KappaType, Field>(kappa, phi) {
    // add default boundary handlers for Corrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <typename KappaType, typename NonOrthoCorrector, typename Field>
void Corrected<KappaType, NonOrthoCorrector, Field>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::Corrected", *this);
}

// NOTE: If KappaType is a field::Tensor, then kappa of each face should be transposed before it
// gets multiplied by face.areaVector(). For now, all the tensor kappas we use are diagonal so it
// does not matter.
template <typename KappaType, typename NonOrthoCorrector, typename Field>
void Corrected<KappaType, NonOrthoCorrector, Field>::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field().mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const double d_CF_norm = d_CF.norm();
    const Vector3d e = d_CF / d_CF_norm;

    const auto Sf = mesh::outwardAreaVector(face, owner);
    const Vector3d Sf_prime = detail::valueAtFace(this->kappa(), face) * Sf;
    const auto [Ef_prime, Tf_prime] = _corrector.decompose(Sf_prime, e);

    // geometric diffusion coefficient
    const double g_diff = Ef_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

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
    const Vector3d grad_f = this->field().gradAtFace(face);
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}

} // namespace prism::scheme::diffusion
