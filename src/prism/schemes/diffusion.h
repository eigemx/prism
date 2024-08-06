#pragma once

#include <cassert>

#include "boundary.h"
#include "diffusion_boundary.h"
#include "fvscheme.h"
#include "prism/field/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/nonortho/nonortho.h"
#include "prism/types.h"

namespace prism::scheme::diffusion {

// Abstract base for diffusion schemes
class IDiffusion {};

template <typename KappaType = field::UniformScalar,
          typename NonOrthoCorrector = nonortho::OverRelaxedCorrector,
          typename GradScheme = gradient::LeastSquares,
          typename Field = field::Scalar>
class CorrectedDiffusion : public IDiffusion, public FVScheme<Field> {
  public:
    CorrectedDiffusion(KappaType kappa, Field phi);

    void apply() override;
    auto field() -> std::optional<Field> override { return _phi; }
    auto corrector() -> NonOrthoCorrector { return _corrector; }
    auto grad_scheme() -> GradScheme& { return _grad_scheme; }
    auto kappa() -> KappaType { return _kappa; }

    using Scheme = CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>;
    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<Scheme,
                                                 boundary::FVSchemeBoundaryHandler<Scheme>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    // TODO: remove this
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    Field _phi;
    KappaType _kappa;
    NonOrthoCorrector _corrector;
    GradScheme _grad_scheme;
    BoundaryHandlersManager _bc_manager;
};

template <typename KappaType = field::UniformScalar, typename Field = field::Scalar>
class NonCorrectedDiffusion : public IDiffusion, public FVScheme<Field> {
  public:
    NonCorrectedDiffusion(KappaType kappa, field::Scalar phi);

    void apply() override;
    auto field() -> std::optional<field::Scalar> override { return _phi; }
    auto needsCorrection() const -> bool override { return false; }
    auto kappa() -> KappaType { return _kappa; }

    using Scheme = NonCorrectedDiffusion<KappaType, Field>;
    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<Scheme,
                                                 boundary::FVSchemeBoundaryHandler<Scheme>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    Field _phi;
    KappaType _kappa;
    BoundaryHandlersManager _bc_manager;
};

//
// CorrectedDiffusion implementation
//
template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>::CorrectedDiffusion(
    KappaType kappa,
    Field phi)
    : _phi(phi), FVScheme<Field>(phi.mesh().nCells()), _kappa(kappa), _grad_scheme(phi) {
    assert(this->needsCorrection() == true &&
           "CorrectedDiffusion::needsCorrection() must return true");

    // add default boundary handlers for CorrectedDiffusion
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Symmetry<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::FixedGradient<Scheme>>();
}

template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
void inline CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied per face basis, using apply_interior() and
     * apply_boundary() functions.
     *
     * The function will check at first if the scheme has completed the first iteration,
     * and if the scheme does not require explicit correction, it will not re-calculate
     * the scheme coefficients and will not zero out the scheme matrix and RHS vector.
     */

    apply_boundary();

    for (const auto& iface : _phi.mesh().interiorFaces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
void inline CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>::
    apply_boundary() {
    boundary::detail::apply_boundary("CorrectedDiffusion", *this);
}

template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
void inline CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>::apply_interior(
    const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const mesh::Cell& neighbor = _phi.mesh().cell(face.neighbor().value());

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const double d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    const Vector3d e = d_CF / d_CF_norm;
    const Vector3d Sf_prime = _kappa.valueAtFace(face) * face.area_vector();

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
    const Vector3d grad_f = _grad_scheme.gradient_at_face(face);

    // update right hand side
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}


//
// NonCorrectedDiffusion implementation
//
template <typename KappaType, typename Field>
NonCorrectedDiffusion<KappaType, Field>::NonCorrectedDiffusion(KappaType kappa, field::Scalar phi)
    : _phi(phi), FVScheme<Field>(phi.mesh().nCells()), _kappa(kappa) {
    assert(this->needsCorrection() == false &&
           "NonCorrectedDiffusion::needsCorrection() must return false");

    // add default boundary handlers for NonCorrectedDiffusion
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Symmetry<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::FixedGradient<Scheme>>();
}

template <typename KappaType, typename Field>
void inline NonCorrectedDiffusion<KappaType, Field>::apply() {
    apply_boundary();

    for (const auto& iface : _phi.mesh().interiorFaces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

template <typename KappaType, typename Field>
void inline NonCorrectedDiffusion<KappaType, Field>::apply_boundary() {
    boundary::detail::apply_boundary("NonCorrectedDiffusion", *this);
}

template <typename KappaType, typename Field>
void inline NonCorrectedDiffusion<KappaType, Field>::apply_interior(const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const mesh::Cell& neighbor = _phi.mesh().cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.area_vector();
    Vector3d Sf_prime = _kappa.valueAtFace(face) * Sf;

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
