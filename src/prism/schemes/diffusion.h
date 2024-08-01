#pragma once

#include <cassert>
#include <memory>
#include <type_traits>

#include "boundary.h"
#include "diffusion_boundary.h"
#include "fmt/core.h"
#include "fvscheme.h"
#include "prism/exceptions.h"
#include "prism/field/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/nonortho/nonortho.h"
#include "prism/types.h"
#include "spdlog/spdlog.h"

namespace prism::diffusion {

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

    using BCManager = boundary::BoundaryHandlersManager<
        CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>>;
    auto bc_manager() -> BCManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    // TODO: remove this
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    const Field _phi;
    const KappaType _kappa;
    NonOrthoCorrector _corrector;
    GradScheme _grad_scheme;
    BCManager _bc_manager;
};

template <typename KappaType = field::UniformScalar, typename Field = field::Scalar>
class NonCorrectedDiffusion : public IDiffusion, public FVScheme<Field> {
  public:
    NonCorrectedDiffusion(KappaType kappa, field::Scalar phi);

    void apply() override;
    auto field() -> std::optional<field::Scalar> override { return _phi; }
    auto requires_correction() const -> bool override { return false; }
    auto kappa() -> KappaType { return _kappa; }

    using BCManager = boundary::BoundaryHandlersManager<NonCorrectedDiffusion<KappaType, Field>>;
    auto bc_manager() -> BCManager& { return _bc_manager; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    const Field _phi;
    const KappaType _kappa;
    BCManager _bc_manager;
};

//
// CorrectedDiffusion implementation
//
template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
CorrectedDiffusion<KappaType, NonOrthoCorrector, GradScheme, Field>::CorrectedDiffusion(
    KappaType kappa,
    Field phi)
    : _phi(phi), FVScheme<Field>(phi.mesh().n_cells()), _kappa(kappa), _grad_scheme(phi) {
    assert(this->requires_correction() == true &&
           "CorrectedDiffusion::requires_correction() must return true");

    // add default boundary handlers for CorrectedDiffusion
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<boundary::Symmetry<Scheme>>();
    _bc_manager.template add_handler<boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<boundary::FixedGradient<Scheme>>();
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

    for (const auto& iface : _phi.mesh().interior_faces()) {
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

    const auto& [Sf, Ef, Tf] = _corrector.interior_triplet(owner, neighbor, face);
    Vector3d Ef_prime = _kappa.value_at_face(face) * Ef;

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
    Vector3d Tf_prime = _kappa.value_at_face(face) * Tf;

    // update right hand side
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}


//
// NonCorrectedDiffusion implementation
//
template <typename KappaType, typename Field>
NonCorrectedDiffusion<KappaType, Field>::NonCorrectedDiffusion(KappaType kappa, field::Scalar phi)
    : _phi(phi), FVScheme<Field>(phi.mesh().n_cells()), _kappa(kappa) {
    assert(this->requires_correction() == false &&
           "NonCorrectedDiffusion::requires_correction() must return false");

    // add default boundary handlers for NonCorrectedDiffusion
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<boundary::Symmetry<Scheme>>();
    _bc_manager.template add_handler<boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<boundary::FixedGradient<Scheme>>();
}

template <typename KappaType, typename Field>
void inline NonCorrectedDiffusion<KappaType, Field>::apply() {
    apply_boundary();

    for (const auto& iface : _phi.mesh().interior_faces()) {
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
    Vector3d Sf_prime = _kappa.value_at_face(face) * Sf;

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
} // namespace prism::diffusion
