#pragma once

#include <cassert>
#include <memory>
#include <type_traits>

#include "boundary.h"
#include "fmt/core.h"
#include "fvscheme.h"
#include "prism/exceptions.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/nonortho/nonortho.h"
#include "prism/types.h"
#include "spdlog/spdlog.h"

namespace prism::diffusion {
//
// forward declarations
//

// Abstract base for diffusion schemes
class IDiffusion;

template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
class CorrectedDiffusion;

template <typename KappaType, typename Field>
class NonCorrectedDiffusion;
} // namespace prism::diffusion

namespace prism::boundary {

//
// CorrectedDiffusion default boundary handlers
//

// Boundary handler for Fixed bounndary condition defined for CorrectedDiffusion with a conserved
// general scalar field (for example: velocity component or temperature).
template <typename K, typename N, typename G>
class Fixed<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>> {
  public:
    void apply(diffusion::CorrectedDiffusion<K, N, G, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed"; }
};

// Boundary handler for symmetry boundary condition, or zero gradient boundary condition. This is
// a special case of the general Neumann boundary condition, where the gradient of the field is
// zero at the boundary (flux is zero), and will not result in any contribution to the right hand
// side of the equation, or the matrix coefficients, and no need for non-orthogonal correction.
// check equation 8.41 - Chapter 8 (Moukallad et al., 2015) and the following paragraph, and
// paragraph 8.6.8.2 - Chapter 8 in same reference.
template <>
class Symmetry<diffusion::IDiffusion> : public FVSchemeBoundaryHandler<diffusion::IDiffusion> {
    void apply(diffusion::IDiffusion& scheme, const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

// We treat outlet boundary condition in diffusion scheme same as symmetry (zero gradient of the
// conserved scalar field)
template <>
class Outlet<diffusion::IDiffusion> : public FVSchemeBoundaryHandler<diffusion::IDiffusion> {
    void apply(diffusion::IDiffusion& scheme, const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "outlet"; }
};

// general Von Neumann boundary condition, or fixed gradient boundary condition.
template <typename K, typename N, typename G>
class FixedGradient<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>> {
  public:
    void apply(diffusion::CorrectedDiffusion<K, N, G, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

//
// NonCorrectedDiffusion default boundary handlers
//
template <typename K>
class Fixed<diffusion::NonCorrectedDiffusion<K, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::NonCorrectedDiffusion<K, field::Scalar>> {
  public:
    void apply(diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename K>
class Symmetry<diffusion::NonCorrectedDiffusion<K, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::NonCorrectedDiffusion<K, field::Scalar>> {
  public:
    void apply(diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename K>
class Outlet<diffusion::NonCorrectedDiffusion<K, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::NonCorrectedDiffusion<K, field::Scalar>> {
  public:
    void apply(diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override {}
    auto inline name() const -> std::string override { return "outlet"; }
};

template <typename K>
class FixedGradient<diffusion::NonCorrectedDiffusion<K, field::Scalar>>
    : public FVSchemeBoundaryHandler<diffusion::NonCorrectedDiffusion<K, field::Scalar>> {
  public:
    void apply(diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

} // namespace prism::boundary

namespace prism::diffusion {

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

    using Scheme = std::remove_reference_t<decltype(*this)>;

    // clang-format off
    _bc_manager.add_handler("empty", 
                            &boundary::create_handler_instance<boundary::Empty<Scheme>>);
    _bc_manager.add_handler("fixed", 
                            &boundary::create_handler_instance<boundary::Fixed<Scheme>>);
    _bc_manager.add_handler("symmetry",
                            &boundary::create_handler_instance<boundary::Symmetry<Scheme>>);
    _bc_manager.add_handler("outlet",
                            &boundary::create_handler_instance<boundary::Outlet<Scheme>>);
    _bc_manager.add_handler("fixed-gradient",
                            &boundary::create_handler_instance<boundary::FixedGradient<Scheme>>);
    // clang-format on
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
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundary_patches()) {
        const mesh::BoundaryCondition& bc = patch.get_bc(_phi.name());
        spdlog::debug(
            "CorrectedDiffusion::apply_boundary(): assigning a boundary handler for boundary "
            "condition type '{}' in patch '{}'.",
            bc.kind_string(),
            patch.name());

        auto handler_creator_opt = _bc_manager.get_handler(bc.kind_string());

        if (!handler_creator_opt.has_value()) {
            throw error::NonImplementedBoundaryCondition(
                "CorrectedDiffusion::apply_boundary()", patch.name(), bc.kind_string());
        }

        auto handler = handler_creator_opt.value()();

        using Scheme = std::remove_reference_t<decltype(*this)>;
        auto fv_handler =
            std::dynamic_pointer_cast<boundary::FVSchemeBoundaryHandler<Scheme>>(handler);

        spdlog::debug(
            "CorrectedDiffusion::apply_boundary(): Applying boundary condition type '{}' on "
            "patch '{}'.",
            fv_handler->name(),
            patch.name());

        fv_handler->apply(*this, patch);
    }
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

    // kappa * g_diff * (Φ_C - Φ_N)
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

    using Scheme = std::remove_reference_t<decltype(*this)>;

    // clang-format off
    _bc_manager.add_handler("empty", 
                            &boundary::create_handler_instance<boundary::Empty<Scheme>>);
    _bc_manager.add_handler("fixed", 
                            &boundary::create_handler_instance<boundary::Fixed<Scheme>>);
    _bc_manager.add_handler("symmetry",
                            &boundary::create_handler_instance<boundary::Symmetry<Scheme>>);
    _bc_manager.add_handler("outlet",
                            &boundary::create_handler_instance<boundary::Outlet<Scheme>>);
    _bc_manager.add_handler("fixed-gradient",
                            &boundary::create_handler_instance<boundary::FixedGradient<Scheme>>);
    // clang-format on
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
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundary_patches()) {
        const mesh::BoundaryCondition& bc = patch.get_bc(_phi.name());
        spdlog::debug(
            "NonCorrectedDiffusion::apply_boundary(): assigning a boundary handler for boundary "
            "condition type '{}' in patch '{}'.",
            bc.kind_string(),
            patch.name());

        auto handler_creator_opt = _bc_manager.get_handler(bc.kind_string());

        if (!handler_creator_opt.has_value()) {
            throw error::NonImplementedBoundaryCondition(
                "NonCorrectedDiffusion::apply_boundary()", patch.name(), bc.kind_string());
        }

        auto handler = handler_creator_opt.value()();

        using Scheme = std::remove_reference_t<decltype(*this)>;
        auto fv_handler =
            std::dynamic_pointer_cast<boundary::FVSchemeBoundaryHandler<Scheme>>(handler);

        spdlog::debug(
            "NonCorrectedDiffusion::apply_boundary(): Applying boundary condition type '{}' on "
            "patch '{}'.",
            fv_handler->name(),
            patch.name());

        fv_handler->apply(*this, patch);
    }
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

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);
}
} // namespace prism::diffusion

namespace prism::boundary {
template <typename K, typename N, typename G>
void Fixed<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>>::apply(
    diffusion::CorrectedDiffusion<K, N, G, field::Scalar>& scheme,
    const mesh::BoundaryPatch& patch) const {
    assert(scheme.field().has_value());
    const auto phi = scheme.field().value();
    const auto& mesh = phi.mesh();
    const auto& corrector = scheme.corrector();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.faces_ids()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        // get the fixed phi variable associated with the face
        const double phi_wall = phi.value_at_face(face);

        const std::size_t cell_id = owner.id();

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();
        const Vector3d e = d_Cf / d_Cf_norm;

        const auto& [_, Ef, Tf] = corrector.boundary_triplet(owner, face);
        Vector3d Ef_prime = kappa.value_at_cell(owner) * Ef;

        const double g_diff = Ef_prime.norm() / (d_Cf_norm + EPSILON);

        scheme.insert(cell_id, cell_id, g_diff);
        scheme.rhs(cell_id) += g_diff * phi_wall;

        Vector3d Tf_prime = kappa.value_at_cell(owner) * Tf;

        // correct non-orhtogonality
        const double phi_c = phi.value_at_cell(owner);
        auto grad_f = ((phi_wall - phi_c) / (d_Cf_norm + EPSILON)) * e;
        scheme.rhs(owner.id()) += Tf_prime.dot(grad_f);
    }
}

template <typename K, typename N, typename G>
void FixedGradient<diffusion::CorrectedDiffusion<K, N, G, field::Scalar>>::apply(
    diffusion::CorrectedDiffusion<K, N, G, field::Scalar>& scheme,
    const mesh::BoundaryPatch& patch) const {
    /** @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face, and the boundary condition
     * is a general Von Neumann boundary condition, or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    const auto phi = scheme.field().value();
    const auto& kappa = scheme.kappa();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.faces_ids()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());

        // get the fixed gradient (flux) value associated with the face
        const auto& boundary_patch = mesh.face_boundary_patch(face);
        const Vector3d wall_grad = boundary_patch.get_vector_bc(phi.name());

        const Vector3d& Sf = face.area_vector();
        Vector3d Sf_prime = kappa.value_at_cell(owner) * Sf;

        // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
        // and paragraph 8.6.8.2
        scheme.rhs(owner.id()) += wall_grad.dot(Sf_prime);
    }
}

template <typename K>
void Fixed<diffusion::NonCorrectedDiffusion<K, field::Scalar>>::apply(
    diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
    const mesh::BoundaryPatch& patch) const {
    assert(scheme.field().has_value());
    const auto phi = scheme.field().value();
    const auto& mesh = phi.mesh();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.faces_ids()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        // get the fixed phi variable associated with the face
        const double phi_wall = phi.value_at_face(face);

        const std::size_t cell_id = owner.id();

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();
        const Vector3d e = d_Cf / d_Cf_norm;

        Vector3d Sf_prime = kappa.value_at_cell(owner) * face.area_vector();

        const double g_diff = Sf_prime.norm() / (d_Cf_norm + EPSILON);

        scheme.insert(cell_id, cell_id, g_diff);
        scheme.rhs(cell_id) += g_diff * phi_wall;
    }
}

template <typename K>
void FixedGradient<diffusion::NonCorrectedDiffusion<K, field::Scalar>>::apply(
    diffusion::NonCorrectedDiffusion<K, field::Scalar>& scheme,
    const mesh::BoundaryPatch& patch) const {
    // This is exactly the same implementation of FixedGradient for CorrectedDiffusion.
    const auto phi = scheme.field().value();
    const auto& kappa = scheme.kappa();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.faces_ids()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());

        // get the fixed gradient (flux) value associated with the face
        const auto& boundary_patch = mesh.face_boundary_patch(face);
        const Vector3d wall_grad = boundary_patch.get_vector_bc(phi.name());

        const Vector3d& Sf = face.area_vector();
        Vector3d Sf_prime = kappa.value_at_cell(owner) * Sf;

        // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
        // and paragraph 8.6.8.2
        scheme.rhs(owner.id()) += wall_grad.dot(Sf_prime);
    }
}

} // namespace prism::boundary
