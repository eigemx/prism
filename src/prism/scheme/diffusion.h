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

// Simple base type for all diffusion schemes, without templating clutter.
class IDiffusion {};

// Base class for all diffusion schemes with a predefined apply() method.
template <typename Kappa, typename Field>
class IAppliedDiffusion : public IFullScheme<Field> {
  public:
    IAppliedDiffusion(Kappa kappa, Field phi);
    void apply() override;
    auto field() -> Field override { return _phi; }
    auto kappa() -> Kappa { return _kappa; }

    using FieldType = Field;
    using KappaType = Kappa;

  private:
    virtual void applyInterior(const mesh::Face& face) = 0;
    virtual void applyBoundary() = 0;
    Field _phi;
    Kappa _kappa;
};

// Base class for all diffusion schemes that need non-orthogonal correction.
class ICorrected : public IDiffusion {};

// Base class for all diffusion schemes that do not need non-orthogonal correction.
class INonCorrected : public IDiffusion {};

// Implementation of diffusion scheme with non-orthogonal correction.
template <typename KappaType = field::UniformScalar,
          typename NonOrthoCorrector = nonortho::OverRelaxedCorrector,
          typename Field = field::Scalar>
class Corrected
    : public ICorrected,
      public IAppliedDiffusion<KappaType, Field>,
      public prism::boundary::BHManagersProvider<
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

// Implementation of diffusion scheme without non-orthogonal correction.
template <typename KappaType = field::UniformScalar, typename Field = field::Scalar>
class NonCorrected : public INonCorrected,
                     public IAppliedDiffusion<KappaType, Field>,
                     public prism::boundary::BHManagersProvider<
                         boundary::ISchemeBoundaryHandler<NonCorrected<KappaType, Field>>> {
  public:
    NonCorrected(KappaType kappa, Field phi);
    auto needsCorrection() const noexcept -> bool override { return false; }

    using Scheme = NonCorrected<KappaType, Field>;

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

//
// diffusion::IAppliedDiffusion implementation
//
template <typename KappaType, typename Field>
IAppliedDiffusion<KappaType, Field>::IAppliedDiffusion(KappaType kappa, Field phi)
    : _kappa(kappa), _phi(phi), IFullScheme<Field>(phi.mesh().cellCount()) {}

template <typename KappaType, typename Field>
void IAppliedDiffusion<KappaType, Field>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied using applyInterior() (per face basis) and
     * applyBoundary() functions.
     *
     */
    applyBoundary();

    const auto& interior_faces = this->field().mesh().interiorFaces();
    std::for_each(interior_faces.begin(), interior_faces.end(), [this](const mesh::Face& face) {
        applyInterior(face);
    });

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

//
// diffusion::Corrected implementation
//
template <typename KappaType, typename NonOrthoCorrector, typename Field>
Corrected<KappaType, NonOrthoCorrector, Field>::Corrected(KappaType kappa, Field phi)
    : IAppliedDiffusion<KappaType, Field>(kappa, phi) {
    // add default boundary handlers for Corrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Empty<Scheme>>();
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
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // vector joining the centers of the two cells
    const Vector3d dCF = neighbor.center() - owner.center();
    const double dCF_norm = dCF.norm();
    const Vector3d e = dCF / dCF_norm;

    const auto Sf = mesh::outwardAreaVector(face, owner);
    const Vector3d Sf_prime = this->kappa().valueAtFace(face) * Sf;
    const auto [Ef_prime, Tf_prime] = _corrector.decompose(Sf_prime, e);

    // geometric diffusion coefficient
    const double gdiff = Ef_prime.norm() / (dCF_norm + EPSILON);

    // gdiff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, gdiff);
    this->insert(neighbor_id, neighbor_id, gdiff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -gdiff);
    this->insert(neighbor_id, owner_id, -gdiff);

    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukallad et al., 2015)
    const Vector3d grad_f = this->field().gradAtFace(face);

    // update right hand side
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}

//
// diffusion::NonCorrected implementation
//
template <typename KappaType, typename Field>
NonCorrected<KappaType, Field>::NonCorrected(KappaType kappa, Field phi)
    : IAppliedDiffusion<KappaType, Field>(kappa, phi) {
    assert(this->needsCorrection() == false &&
           "prism::scheme::diffusion::NonCorrected::needsCorrection() must return false");

    // add default boundary handlers for NonCorrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Empty<Scheme>>();
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
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.areaVector();
    Vector3d Sf_prime = this->kappa().valueAtFace(face) * Sf;

    // geometric diffusion coefficient
    const double g_diff = Sf_prime.norm() / (d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);
}
} // namespace prism::scheme::diffusion
