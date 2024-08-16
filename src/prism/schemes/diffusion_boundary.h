#pragma once

#include "boundary.h"
#include "prism/constants.h"
#include "prism/field/scalar.h"

namespace prism::scheme::diffusion {

class ICorrected;

template <typename T>
concept ICorrectedBased = std::is_base_of_v<ICorrected, T>;

class INonCorrected;

template <typename T>
concept INonCorrectedBased = std::is_base_of_v<INonCorrected, T>;

template <typename KappaType, typename NonOrthoCorrector, typename GradScheme, typename Field>
class Corrected;

template <typename KappaType, typename Field>
class NonCorrected;


} // namespace prism::scheme::diffusion

namespace prism::scheme::boundary {
//
// diffusion::Corrected default boundary handlers
//
// Boundary handler for symmetry boundary condition, or zero gradient boundary condition. This is
// a special case of the general Neumann boundary condition, where the gradient of the field is
// zero at the boundary (flux is zero), and will not result in any contribution to the right hand
// side of the equation, or the matrix coefficients, and no need for non-orthogonal correction.
// check equation 8.41 - Chapter 8 (Moukallad et al., 2015) and the following paragraph, and
// paragraph 8.6.8.2 - Chapter 8 in same reference.
template <scheme::diffusion::ICorrectedBased Scheme>
class Symmetry<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

// We treat outlet boundary condition in diffusion scheme same as symmetry (zero gradient of the
// conserved scalar field)
template <scheme::diffusion::ICorrectedBased Scheme>
class Outlet<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "outlet"; }
};

// Boundary handler for Fixed bounndary condition defined for CorrectedDiffusion with a conserved
// general scalar field (for example: velocity component or temperature).
template <scheme::diffusion::ICorrectedBased Scheme>
class Fixed<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <scheme::diffusion::ICorrectedBased Scheme>
class NoSlip<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

// general Von Neumann boundary condition, or fixed gradient boundary condition.
template <scheme::diffusion::ICorrectedBased Scheme>
class FixedGradient<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

//
// diffusion::NonCorrected default boundary handlers
//
template <scheme::diffusion::INonCorrectedBased Scheme>
class Fixed<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <scheme::diffusion::INonCorrectedBased Scheme>
class Symmetry<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <scheme::diffusion::INonCorrectedBased Scheme>
class Outlet<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "outlet"; }
};

template <scheme::diffusion::INonCorrectedBased Scheme>
class FixedGradient<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

template <scheme::diffusion::INonCorrectedBased Scheme>
class NoSlip<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

// TODO: boundary handlers for Corrected and NonCorrected should be the same,
// right?
template <scheme::diffusion::ICorrectedBased Scheme>
void Fixed<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    assert(scheme.field().has_value());
    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();
    const auto& corrector = scheme.corrector();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        // get the fixed phi variable associated with the face
        const double phi_wall = phi.valueAtFace(face);

        const std::size_t cell_id = owner.id();

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();
        const Vector3d e = d_Cf / d_Cf_norm;

        const Vector3d Sf_prime = kappa.valueAtFace(face) * face.area_vector();

        const auto& [Ef_prime, Tf_prime] = corrector.decompose(Sf_prime, e);

        const double g_diff = Ef_prime.norm() / (d_Cf_norm + EPSILON);

        scheme.insert(cell_id, cell_id, g_diff);
        scheme.rhs(cell_id) += g_diff * phi_wall;

        // correct non-orhtogonality
        const double phi_c = phi.valueAtCell(owner);
        auto grad_f = ((phi_wall - phi_c) / (d_Cf_norm + EPSILON)) * e;
        scheme.rhs(owner.id()) += Tf_prime.dot(grad_f);
    }
}

template <scheme::diffusion::ICorrectedBased Scheme>
void NoSlip<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    Fixed<Scheme> fixed;
    return fixed.apply(scheme, patch);
}

template <scheme::diffusion::ICorrectedBased Scheme>
void FixedGradient<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face, and the boundary condition
     * is a general Von Neumann boundary condition, or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    const auto phi = scheme.field();
    const auto& kappa = scheme.kappa();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());

        // get the fixed gradient (flux) value associated with the face
        const auto& boundary_patch = mesh.faceBoundaryPatch(face);
        const Vector3d wall_grad = boundary_patch.getVectorBoundaryCondition(phi.name());

        const Vector3d& Sf = face.area_vector();
        Vector3d Sf_prime = kappa.valueAtCell(owner) * Sf;

        // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
        // and paragraph 8.6.8.2
        scheme.rhs(owner.id()) += wall_grad.dot(Sf_prime);
    }
}

template <scheme::diffusion::INonCorrectedBased Scheme>
void Fixed<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    assert(scheme.field().has_value());
    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        // get the fixed phi variable associated with the face
        const double phi_wall = phi.valueAtFace(face);

        const std::size_t cell_id = owner.id();

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();
        const Vector3d e = d_Cf / d_Cf_norm;

        Vector3d Sf_prime = kappa.valueAtCell(owner) * face.area_vector();

        const double g_diff = Sf_prime.norm() / (d_Cf_norm + EPSILON);

        scheme.insert(cell_id, cell_id, g_diff);
        scheme.rhs(cell_id) += g_diff * phi_wall;
    }
}

template <scheme::diffusion::INonCorrectedBased Scheme>
void NoSlip<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    Fixed<Scheme> fixed;
    return fixed.apply(scheme, patch);
}

template <scheme::diffusion::INonCorrectedBased Scheme>
void FixedGradient<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    // This is exactly the same implementation of FixedGradient for CorrectedDiffusion.
    const auto phi = scheme.field();
    const auto& kappa = scheme.kappa();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());

        // get the fixed gradient (flux) value associated with the face
        const auto& boundary_patch = mesh.faceBoundaryPatch(face);
        const Vector3d wall_grad = boundary_patch.getVectorBoundaryCondition(phi.name());

        const Vector3d& Sf = face.area_vector();
        Vector3d Sf_prime = kappa.valueAtCell(owner) * Sf;

        // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
        // and paragraph 8.6.8.2
        scheme.rhs(owner.id()) += wall_grad.dot(Sf_prime);
    }
}

} // namespace prism::scheme::boundary