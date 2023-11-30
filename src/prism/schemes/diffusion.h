#pragma once

#include "fvscheme.h"
#include "prism/exceptions.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/nonortho/nonortho.h"
#include "prism/types.h"

namespace prism::diffusion {

template <typename NonOrthoCorrector = nonortho::OverRelaxedCorrector<>>
class AbstractDiffusion : public FVScheme {
  public:
    AbstractDiffusion(ScalarField phi);

    void apply() override;
    auto field() -> std::optional<ScalarField> override { return _phi; }

    // This will be overloaded to return false in case of NilCorrector
    auto requires_correction() const -> bool override { return true; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override;

    virtual auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Face& face)
        -> Vector3d = 0;
    virtual auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Cell& cell)
        -> Vector3d = 0;

    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void apply_boundary_gradient(const mesh::Cell& cell, const mesh::Face& face);
    void correct_nonorhto_boundary_fixed(const mesh::Cell& cell,
                                         const mesh::Face& face,
                                         const Vector3d& Tf_prime);

    const ScalarField _phi;
    NonOrthoCorrector _corrector;
};

template <typename DiffusionCoeffType, typename Corrector>
class Diffusion : public AbstractDiffusion<Corrector> {
  public:
    Diffusion(DiffusionCoeffType kappa, ScalarField phi);

  private:
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Face& face)
        -> Vector3d override;
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Cell& cell)
        -> Vector3d override;
    DiffusionCoeffType _kappa;
};

template <typename Corrector>
class Diffusion<double, Corrector> : public AbstractDiffusion<Corrector> {
  public:
    Diffusion(double kappa, ScalarField phi);

  private:
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Face& face)
        -> Vector3d override;
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Cell& cell)
        -> Vector3d override;
    double _kappa;
};

template <typename Corrector>
class Diffusion<Matrix3d, Corrector> : public AbstractDiffusion<Corrector> {
  public:
    Diffusion(Matrix3d kappa, ScalarField phi);

  private:
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Face& face)
        -> Vector3d override;
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Cell& cell)
        -> Vector3d override;
    Matrix3d _kappa;
};

template <typename Corrector>
class Diffusion<TensorField, Corrector> : public AbstractDiffusion<Corrector> {
  public:
    Diffusion(TensorField kappa, ScalarField phi);

  private:
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Face& face)
        -> Vector3d override;
    auto diffusion_scaled_vector(const Vector3d& vec, const mesh::Cell& cell)
        -> Vector3d override;
    TensorField _kappa;
};

template <typename Corrector>
Diffusion<double, Corrector>::Diffusion(double kappa, ScalarField phi)
    : _kappa(kappa), AbstractDiffusion<Corrector>(std::move(phi)) {}

template <typename Corrector>
Diffusion<Matrix3d, Corrector>::Diffusion(Matrix3d kappa, ScalarField phi)
    : _kappa(std::move(kappa)), AbstractDiffusion<Corrector>(std::move(phi)) {}

template <typename Corrector>
Diffusion<TensorField, Corrector>::Diffusion(TensorField kappa, ScalarField phi)
    : _kappa(std::move(kappa)), AbstractDiffusion<Corrector>(std::move(phi)) {}

template <typename NonOrthoCorrector>
AbstractDiffusion<NonOrthoCorrector>::AbstractDiffusion(ScalarField phi)
    : _phi(phi), FVScheme(phi.mesh().n_cells()), _corrector(phi) {}

template <>
AbstractDiffusion<nonortho::NilCorrector>::AbstractDiffusion(ScalarField phi) // NOLINT
    : _phi(phi), FVScheme(phi.mesh().n_cells()) {}

template <typename Corrector>
auto Diffusion<double, Corrector>::Diffusion::diffusion_scaled_vector(
    const Vector3d& vec,
    const mesh::Face& face) // NOLINT
    -> Vector3d {
    return _kappa * vec;
}

template <typename Corrector>
auto Diffusion<double, Corrector>::Diffusion::diffusion_scaled_vector(
    const Vector3d& vec,
    const mesh::Cell& cell) // NOLINT
    -> Vector3d {
    return _kappa * vec;
}

template <typename Corrector>
auto Diffusion<Matrix3d, Corrector>::Diffusion::diffusion_scaled_vector(
    const Vector3d& vec,
    const mesh::Face& face) // NOLINT
    -> Vector3d {
    return _kappa * vec;
}

template <typename Corrector>
auto Diffusion<Matrix3d, Corrector>::Diffusion::diffusion_scaled_vector(
    const Vector3d& vec,
    const mesh::Cell& cell) // NOLINT
    -> Vector3d {
    return _kappa * vec;
}

template <typename Corrector>
auto Diffusion<TensorField, Corrector>::Diffusion::diffusion_scaled_vector(const Vector3d& vec,
                                                                           const mesh::Face& face)
    -> Vector3d {
    return _kappa.value_at_face(face) * vec;
}

template <typename Corrector>
auto Diffusion<TensorField, Corrector>::Diffusion::diffusion_scaled_vector(
    const Vector3d& vec,
    const mesh::Cell& cell) // NOLINT
    -> Vector3d {
    return _kappa.value_at_cell(cell) * vec;
}

template <typename NonOrthoCorrector>
void inline AbstractDiffusion<NonOrthoCorrector>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied per face basis, using apply_interior() and 
     * apply_boundary() functions.
     * 
     * The function will check at first if the scheme has completed the first iteration,
     * and if the scheme does not require explicit correction, it will not re-calculate
     * the scheme coefficients and will not zero out the scheme matrix and RHS vector.
     */
    for (const auto& bface : _phi.mesh().boundary_faces()) {
        apply_boundary(bface);
    }

    for (const auto& iface : _phi.mesh().interior_faces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    collect();
}

template <typename NonOrthoCorrector>
void AbstractDiffusion<NonOrthoCorrector>::apply_boundary(const mesh::Face& face) {
    /**
     * @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face. The function iteself does not
     * apply the discretized equation, rather it calls the appropriate function
     * based on the boundary type.
     *
     * @param cell The cell that owns the boundary face.
     * @param face The boundary face.
     */
    const auto& owner = _phi.mesh().cell(face.owner());
    const auto& boundary_patch = _phi.mesh().face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.bc_type()) {
        // empty boundary patch, do nothing
        case mesh::BoundaryConditionType::Empty: {
            return;
        }

        // fixed boundary patch, or Dirichlet boundary condition
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            apply_boundary_fixed(owner, face);
            return;
        }

        // Symmetry boundary patch, or zero gradient boundary condition.
        // This is a special case of the general Neumann boundary condition,
        // where the gradient of the field is zero at the boundary (flux is zero),
        // and will not result in any contribution to the right hand side of the equation,
        // or the matrix coefficients, and no need for non-orthogonal correction.
        // check equation 8.41 - Chapter 8 (Moukallad et al., 2015) and the following paragraph,
        // and paragraph 8.6.8.2 - Chapter 8 in same reference.
        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Outlet: {
            return;
        }

        // general Von Neumann boundary condition, or fixed gradient boundary condition.
        case mesh::BoundaryConditionType::FixedGradient: {
            apply_boundary_gradient(owner, face);
            return;
        }

        default:
            throw error::NonImplementedBoundaryCondition(
                "prism::diffusion::Diffusion::apply_boundary()",
                boundary_patch.name(),
                boundary_condition.bc_type_str());
    }
}

template <typename NonOrthoCorrector>
void AbstractDiffusion<NonOrthoCorrector>::apply_boundary_gradient(const mesh::Cell& cell,
                                                                   const mesh::Face& face) {
    /** @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face, and the boundary condition
     * is a general Von Neumann boundary condition, or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    // get the fixed gradient (flux) value associated with the face
    const auto& boundary_patch = _phi.mesh().face_boundary_patch(face);
    const Vector3d wall_grad = boundary_patch.get_vector_bc(_phi.name());

    const Vector3d& Sf = face.area_vector();
    Vector3d Sf_prime = diffusion_scaled_vector(Sf, cell);

    // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
    // and paragraph 8.6.8.2
    rhs(cell.id()) += wall_grad.dot(Sf_prime);
}


template <typename NonOrthoCorrector>
void AbstractDiffusion<NonOrthoCorrector>::apply_interior(const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const mesh::Cell& neighbor = _phi.mesh().cell(face.neighbor().value());

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const double d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    const Vector3d e = d_CF / d_CF_norm;

    const auto& [Sf, Ef, Tf] = _corrector.interior_triplet(owner, neighbor, face);
    Vector3d Ef_prime = diffusion_scaled_vector(Ef, face);

    // geometric diffusion coefficient
    const double g_diff = Ef_prime.norm() / (d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    insert(owner_id, owner_id, g_diff);
    insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    insert(owner_id, neighbor_id, -g_diff);
    insert(neighbor_id, owner_id, -g_diff);

    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukallad et al., 2015)
    const Vector3d grad_f = _corrector.grad_scheme().gradient_at_face(face);
    Vector3d Tf_prime = diffusion_scaled_vector(Tf, face);

    // update right hand side
    rhs(owner_id) += Tf_prime.dot(grad_f);
    rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}

template <>
void inline AbstractDiffusion<nonortho::NilCorrector>::apply_interior(const mesh::Face& face) {
    const mesh::Cell& owner = _phi.mesh().cell(face.owner());
    const mesh::Cell& neighbor = _phi.mesh().cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.area_vector();
    Vector3d Sf_prime = diffusion_scaled_vector(Sf, face);

    // geometric diffusion coefficient
    const double g_diff = Sf_prime.norm() / (d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    insert(owner_id, owner_id, g_diff);
    insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    insert(owner_id, neighbor_id, -g_diff);
    insert(neighbor_id, owner_id, -g_diff);
}


template <typename NonOrthoCorrector>
void AbstractDiffusion<NonOrthoCorrector>::correct_nonorhto_boundary_fixed(
    const mesh::Cell& cell,
    const mesh::Face& face,
    const Vector3d& Tf_prime) {
    // we need to calculate the gradient of phi at the face
    // first let's calculate the vector joining the face center to the cell center
    const Vector3d d_CF = face.center() - cell.center();
    const double d_CF_norm = d_CF.norm();
    const Vector3d e = d_CF / d_CF_norm;

    const double phi_wall = _phi.value_at_face(face);
    const double phi_c = _phi.value_at_cell(cell);

    auto grad_f = ((phi_wall - phi_c) / (d_CF_norm + EPSILON)) * e;
    rhs(cell.id()) += Tf_prime.dot(grad_f);
}


template <typename NonOrthoCorrector>
void AbstractDiffusion<NonOrthoCorrector>::apply_boundary_fixed(const mesh::Cell& cell,
                                                                const mesh::Face& face) {
    // get the fixed phi variable associated with the face
    const double phi_wall = _phi.value_at_face(face);

    const std::size_t cell_id = cell.id();

    // vector joining the centers of the cell and the face
    const Vector3d d_Cf = face.center() - cell.center();
    const double d_Cf_norm = d_Cf.norm();

    const auto& [_, Ef, Tf] = _corrector.boundary_triplet(cell, face);
    Vector3d Ef_prime = diffusion_scaled_vector(Ef, cell);

    const double g_diff = Ef_prime.norm() / (d_Cf_norm + EPSILON);

    insert(cell_id, cell_id, g_diff);
    rhs(cell_id) += g_diff * phi_wall;

    Vector3d Tf_prime = diffusion_scaled_vector(Tf, cell);
    correct_nonorhto_boundary_fixed(cell, face, Tf_prime);
}

template <>
auto inline AbstractDiffusion<nonortho::NilCorrector>::requires_correction() const -> bool {
    return false;
}

} // namespace prism::diffusion