#pragma once

#include <cassert>

#include "boundary.h"
#include "diffusion_boundary.h"
#include "nonortho.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/types.h"
#include "scheme.h"

namespace prism::scheme::diffusion {

template <typename Kappa, typename Field>
class IDiffusion : public IFullScheme<Field> {
  public:
    IDiffusion(Kappa kappa, Field phi);
    auto field() -> Field override { return _phi; }
    auto kappa() -> Kappa { return _kappa; }

    using FieldType = Field;
    using KappaType = Kappa;

  private:
    Field _phi;
    Kappa _kappa;
};

class ICorrected {};
class INonCorrected {};

// TODO: Simplify the BoundaryHandlerManager declaration mess and find a way to use
// BHManagerProvider without type clutter
template <typename KappaType = field::UniformScalar,
          typename NonOrthoCorrector = nonortho::OverRelaxedCorrector,
          typename Field = field::Scalar,
          typename GradScheme = gradient::LeastSquares<Field>>
class Corrected : public ICorrected,
                  public IDiffusion<KappaType, Field>,
                  public gradient::GradientProvider<Field, GradScheme> {
  public:
    Corrected(KappaType kappa, Field phi);

    void apply() override;
    auto corrector() -> NonOrthoCorrector { return _corrector; }
    auto needsCorrection() const noexcept -> bool override { return true; }

    using Scheme = Corrected<KappaType, NonOrthoCorrector, Field, GradScheme>;
    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<boundary::ISchemeBoundaryHandler<Scheme>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    // TODO: remove this
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    NonOrthoCorrector _corrector;
    BoundaryHandlersManager _bc_manager;
};

template <typename KappaType = field::UniformScalar, typename Field = field::Scalar>
class NonCorrected : public INonCorrected, public IDiffusion<KappaType, Field> {
  public:
    NonCorrected(KappaType kappa, Field phi);

    void apply() override;
    auto needsCorrection() const noexcept -> bool override { return false; }

    using Scheme = NonCorrected<KappaType, Field>;
    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<boundary::ISchemeBoundaryHandler<Scheme>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    BoundaryHandlersManager _bc_manager;
};

template <typename KappaType, typename Field>
IDiffusion<KappaType, Field>::IDiffusion(KappaType kappa, Field phi)
    : _kappa(kappa), _phi(phi), IFullScheme<Field>(phi.mesh().nCells()) {}

//
// Corrected implementation
//
template <typename KappaType, typename NonOrthoCorrector, typename Field, typename GradScheme>
Corrected<KappaType, NonOrthoCorrector, Field, GradScheme>::Corrected(KappaType kappa, Field phi)
    : IDiffusion<KappaType, Field>(kappa, phi),
      gradient::GradientProvider<Field, GradScheme>(phi) {
    assert(this->needsCorrection() == true &&
           "prism::scheme::diffusion::Corrected::needsCorrection() must return true");

    // add default boundary handlers for Corrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template addHandler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Symmetry<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::FixedGradient<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::ZeroGradient<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::NoSlip<Scheme>>();
}

template <typename KappaType, typename NonOrthoCorrector, typename Field, typename GradScheme>
void inline Corrected<KappaType, NonOrthoCorrector, Field, GradScheme>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied per face basis, using apply_interior() and
     * apply_boundary() functions.
     *
     * The function will check at first if the scheme has completed the first iteration,
     * and if the scheme does not require explicit correction, it will not re-calculate
     * the scheme coefficients and will not zero out the scheme matrix and RHS vector.
     */

    apply_boundary();

    for (const auto& iface : this->field().mesh().interiorFaces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

template <typename KappaType, typename NonOrthoCorrector, typename Field, typename GradScheme>
void inline Corrected<KappaType, NonOrthoCorrector, Field, GradScheme>::apply_boundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::Corrected", *this);
}

template <typename KappaType, typename NonOrthoCorrector, typename Field, typename GradScheme>
void inline Corrected<KappaType, NonOrthoCorrector, Field, GradScheme>::apply_interior(
    const mesh::Face& face) {
    const auto& mesh = this->field().mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const double d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    const Vector3d e = d_CF / d_CF_norm;
    const Vector3d Sf_prime = this->kappa().valueAtFace(face) * face.area_vector();

    const auto [Ef_prime, Tf_prime] = _corrector.decompose(Sf_prime, e);

    // geometric diffusion coefficient
    const double g_diff = Ef_prime.norm() / (d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);

    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukallad et al., 2015)
    const Vector3d grad_f = this->gradScheme().gradAtFace(face);

    // update right hand side
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}


//
// NonCorrected implementation
//
template <typename KappaType, typename Field>
NonCorrected<KappaType, Field>::NonCorrected(KappaType kappa, Field phi)
    : IDiffusion<KappaType, Field>(kappa, phi) {
    assert(this->needsCorrection() == false &&
           "prism::scheme::diffusion::NonCorrected::needsCorrection() must return false");

    // add default boundary handlers for NonCorrected
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template addHandler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Symmetry<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::FixedGradient<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::ZeroGradient<Scheme>>();
    _bc_manager.template addHandler<scheme::boundary::NoSlip<Scheme>>();
}

template <typename KappaType, typename Field>
void inline NonCorrected<KappaType, Field>::apply() {
    apply_boundary();

    for (const auto& iface : this->field().mesh().interiorFaces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

template <typename KappaType, typename Field>
void inline NonCorrected<KappaType, Field>::apply_boundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::NonCorrected", *this);
}

template <typename KappaType, typename Field>
void inline NonCorrected<KappaType, Field>::apply_interior(const mesh::Face& face) {
    const auto& mesh = this->field().mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.area_vector();
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
