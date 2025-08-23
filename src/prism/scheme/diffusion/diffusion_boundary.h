#pragma once

#include <concepts>

#include "prism/constants.h"
#include "prism/field/ifield.h"
#include "prism/field/tensor.h"
#include "prism/scheme/boundary.h"

namespace prism::scheme::diffusion {

namespace detail {

// The following functions are used to get the value of a field at a face or cell.
// They are specialized for field::Tensor to return a transposed value of the tensor value at the
// face or cell, to follow equation (8.93) from Moukalled et al. (2015). And for other field
// types, they return the value directly.
template <field::IFieldBased Field>
auto valueAtFace(const SharedPtr<Field>& field, const mesh::Face& face) -> Field::ValueType;

template <field::IFieldBased Field>
auto valueAtCell(const SharedPtr<Field>& field, const mesh::Cell& cell) -> Field::ValueType;

template <>
auto valueAtFace(const SharedPtr<field::Tensor>& field, const mesh::Face& face)
    -> field::Tensor::ValueType;

template <>
auto valueAtCell(const SharedPtr<field::Tensor>& field, const mesh::Cell& cell)
    -> field::Tensor::ValueType;
} // namespace detail

// Forward declarations for diffusion scheme classes defined in diffusion.h
class IDiffusion;

template <typename T>
concept IDiffusionBased = std::derived_from<T, IDiffusion>;

class ICorrected;

template <typename T>
concept ICorrectedBased = std::derived_from<T, ICorrected>;

class INonCorrected;

template <typename T>
concept INonCorrectedBased = std::derived_from<T, INonCorrected>;

template <typename Kappa, typename NonOrthoCorrector>
class Corrected;

template <typename Kappa>
class NonCorrected;

} // namespace prism::scheme::diffusion

namespace prism::scheme::boundary {
template <diffusion::IDiffusionBased Scheme>
class Symmetry<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <diffusion::IDiffusionBased Scheme>
class Outlet<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "outlet"; }
};

// general Von Neumann boundary condition, or fixed gradient boundary condition.
template <diffusion::IDiffusionBased Scheme>
class FixedGradient<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

template <diffusion::IDiffusionBased Scheme>
class ZeroGradient<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

//
// diffusion::NonCorrected default boundary handlers
//
template <diffusion::INonCorrectedBased Scheme>
class Fixed<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <diffusion::INonCorrectedBased Scheme>
class NoSlip<Scheme> : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

//
// diffusion::Corrected default boundary handlers
//
// Boundary handler for symmetry boundary condition, or zero gradient boundary condition. This is
// a special case of the general Neumann boundary condition, where the gradient of the field is
// zero at the boundary (flux is zero), and will not result in any contribution to the right hand
// side of the equation, or the matrix coefficients, and no need for non-orthogonal correction.
// check equation 8.41 - Chapter 8 (Moukalled et al., 2015) and the following paragraph, and
// paragraph 8.6.8.2 - Chapter 8 in same reference.

// We treat outlet boundary condition in diffusion scheme same as symmetry (zero gradient of the
// conserved scalar field)


// Boundary handler for Fixed bounndary condition defined for CorrectedDiffusion
template <diffusion::ICorrectedBased Scheme>
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

//
// Implementation of boundary handlers for diffusion schemes
//

template <diffusion::IDiffusionBased Scheme>
void FixedGradient<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation to the cell, when the current face
     * is a boundary face, and the boundary condition is a general Von Neumann boundary condition,
     * or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    const auto& phi = scheme.field();
    const auto& kappa = scheme.kappa();
    const auto& mesh = phi->mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);
        const mesh::Cell& owner = mesh->cell(face.owner());

        // get the fixed gradient (flux) value associated with the face
        const auto& boundary_patch = mesh->faceBoundaryPatch(face);
        const Vector3d wall_grad = boundary_patch.getVectorBoundaryCondition(phi->name());

        const Vector3d& Sf = face.areaVector();
        Vector3d Sf_prime = kappa->valueAtCell(owner) * Sf;

        // check Moukalled et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
        // and paragraph 8.6.8.2
        scheme.rhs(owner.id()) += -wall_grad.dot(Sf_prime);
    }
}

template <scheme::diffusion::INonCorrectedBased Scheme>
void Fixed<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation (without non-orthogonal correction
     * ) to the cell, when the current face is a boundary face, and the boundary condition is a
     * fixed value boundary condition.
     *
     * @param scheme The diffusion scheme.
     * @param patch The boundary patch containing the faces.
     */
    const auto phi = scheme.field();
    const auto& mesh = phi->mesh();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);
        const mesh::Cell& owner = mesh->cell(face.owner());
        const double phi_wall = phi->valueAtFace(face);
        const std::size_t cell_id = owner.id();

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();

        Vector3d Sf_prime = diffusion::detail::valueAtCell(kappa, owner) * face.areaVector();
        const double g_diff = Sf_prime.dot(d_Cf) / (d_Cf_norm * d_Cf_norm + EPSILON);

        scheme.insert(cell_id, cell_id, g_diff);
        scheme.rhs(cell_id) += g_diff * phi_wall;
    }
}

template <scheme::diffusion::INonCorrectedBased Scheme>
void NoSlip<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation to the cell, when the current face
     * is a boundary face, and the boundary condition is a no-slip boundary condition. no-slip
     * boundary condition in diffusion scheme is a special case of the fixed boundary condition.
     *
     * @param scheme The diffusion scheme.
     * @param patch The boundary patch containing the faces.
     */
    Fixed<Scheme> fixed;
    return fixed.apply(scheme, patch);
}

template <scheme::diffusion::ICorrectedBased Scheme>
void Fixed<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation (with non-orthogonal correction)
     * to the cell, when the current face is a boundary face, and the boundary condition is a
     * fixed value boundary condition.
     *
     * @param scheme The diffusion scheme.
     * @param patch The boundary patch containing the faces.
     */
    auto phi = scheme.field();
    const auto& mesh = phi->mesh();
    const auto& corrector = scheme.corrector();
    const auto& kappa = scheme.kappa();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);
        const mesh::Cell& owner = mesh->cell(face.owner());
        const std::size_t owner_id = owner.id();

        // get the fixed phi variable associated with the face
        const double phi_wall = phi->valueAtFace(face);

        // vector joining the centers of the cell and the face
        const Vector3d d_Cf = face.center() - owner.center();
        const double d_Cf_norm = d_Cf.norm();
        const Vector3d e = d_Cf / d_Cf_norm;

        Vector3d Sf_prime = diffusion::detail::valueAtCell(kappa, owner) * face.areaVector();
        const auto& [Ef_prime, Tf_prime] = corrector.decompose(Sf_prime, e);
        const double g_diff = Ef_prime.dot(d_Cf) / (d_Cf_norm * d_Cf_norm + EPSILON);

        scheme.insert(owner_id, owner_id, g_diff);
        scheme.rhs(owner_id) += (g_diff * phi_wall) + Tf_prime.dot(phi->gradAtFace(face));
    }
}

template <scheme::diffusion::ICorrectedBased Scheme>
void NoSlip<Scheme>::apply(Scheme& scheme, const mesh::BoundaryPatch& patch) {
    /** @brief Applies boundary discretized diffusion equation to the cell, when the current face
     * is a boundary face, and the boundary condition is a no-slip boundary condition. no-slip
     * boundary condition in diffusion scheme is a special case of the fixed boundary condition.
     *
     * @param scheme The diffusion scheme.
     * @param patch The boundary patch containing the faces.
     */
    Fixed<Scheme> fixed;
    return fixed.apply(scheme, patch);
}
} // namespace prism::scheme::boundary

namespace prism::scheme::diffusion::detail {
template <field::IFieldBased Field>
auto valueAtFace(const SharedPtr<Field>& field, const mesh::Face& face) -> Field::ValueType {
    return field->valueAtFace(face.id());
}

template <field::IFieldBased Field>
auto valueAtCell(const SharedPtr<Field>& field, const mesh::Cell& cell) -> Field::ValueType {
    return field->valueAtCell(cell.id());
}

template <>
auto inline valueAtFace(const SharedPtr<field::Tensor>& field, const mesh::Face& face)
    -> field::Tensor::ValueType {
    return field->valueAtFace(face.id()).transpose();
}

template <>
auto inline valueAtCell(const SharedPtr<field::Tensor>& field, const mesh::Cell& cell)
    -> field::Tensor::ValueType {
    return field->valueAtCell(cell.id()).transpose();
}
} // namespace prism::scheme::diffusion::detail
