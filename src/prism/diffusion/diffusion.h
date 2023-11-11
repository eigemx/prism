#pragma once

#include "prism/field.h"
#include "prism/fvscheme.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/pmesh.h"
#include "prism/nonortho/nonortho.h"
#include "prism/types.h"

namespace prism::diffusion {

template <typename NonOrthoCorrector = nonortho::OverRelaxedCorrector<>>
class Diffusion : public FVScheme {
  public:
    Diffusion(double kappa, ScalarField& phi);

    void apply() override;
    auto field() -> std::optional<ScalarField> override { return _phi; }
    auto requires_correction() const -> bool override { return true; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void apply_boundary_gradient(const mesh::Cell& cell, const mesh::Face& face);
    void correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                          const mesh::Face& face,
                                          const Vector3d& T_f);

    double _kappa;
    ScalarField _phi;
    const mesh::PMesh& _mesh;
    NonOrthoCorrector _corrector;
};

template <typename NonOrthoCorrector>
Diffusion<NonOrthoCorrector>::Diffusion(double kappa, ScalarField& phi)
    : _kappa(kappa),
      _phi(phi),
      _mesh(phi.mesh()),
      FVScheme(phi.mesh().n_cells()),
      _corrector(phi) {}

template <>
inline Diffusion<nonortho::NilCorrector>::Diffusion(double kappa, ScalarField& phi)
    : _kappa(kappa), _phi(phi), _mesh(phi.mesh()), FVScheme(phi.mesh().n_cells()) {}


template <typename NonOrthoCorrector>
void inline Diffusion<NonOrthoCorrector>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied per face basis, using apply_interior() and 
     * apply_boundary() functions.
     * 
     * The function will check at first if the scheme has completed the first iteration,
     * and if the scheme does not require explicit correction, it will not re-calculate
     * the scheme coefficients and will not zero out the scheme matrix and RHS vector.
     */
    for (const auto& bface : _mesh.boundary_faces()) {
        apply_boundary(_mesh.cell(bface.owner()), bface);
    }

    for (const auto& iface : _mesh.interior_faces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    collect();
}

template <typename NonOrthoCorrector>
void Diffusion<NonOrthoCorrector>::apply_boundary(const mesh::Cell& cell,
                                                  const mesh::Face& face) {
    /**
     * @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face. The function iteself does not
     * apply the discretized equation, rather it calls the appropriate function
     * based on the boundary type.
     *
     * @param cell The cell that owns the boundary face.
     * @param face The boundary face.
     */
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.bc_type()) {
        // empty boundary patch, do nothing
        case mesh::BoundaryConditionType::Empty: {
            return;
        }

        // fixed boundary patch, or Dirichlet boundary condition
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            apply_boundary_fixed(cell, face);
            return;
        }

        // Symmetry boundary patch, or zero gradient boundary condition.
        // This is a special case of the general Neumann boundary condition,
        // where the gradient of the field is zero at the boundary (flux is zero),
        // and will not result in any contribution to the right hand side of the equation,
        // or the matrix coefficients. and no need for non-orthogonal correction.
        // check equation 8.41 - Chapter 8 (Moukallad et al., 2015) and the following paragraph,
        // and paragraph 8.6.8.2 - Chapter 8 in same reference.
        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Outlet: {
            return;
        }

        // general Von Neumann boundary condition, or fixed gradient boundary condition.
        case mesh::BoundaryConditionType::FixedGradient: {
            apply_boundary_gradient(cell, face);
            return;
        }

        default:
            throw std::runtime_error(
                fmt::format("diffusion::Diffusion::apply_boundary(): "
                            "Non-implemented boundary condition type for boundary patch: '{}'",
                            boundary_patch.name()));
    }
}

template <typename NonOrthoCorrector>
void Diffusion<NonOrthoCorrector>::apply_boundary_gradient(const mesh::Cell& cell,
                                                           const mesh::Face& face) {
    /** @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face, and the boundary condition
     * is a general Von Neumann boundary condition, or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    // get the fixed gradient (flux) value associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto flux_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
    // and paragraph 8.6.8.2
    rhs(cell_id) += -flux_wall * face.area();
}


template <typename NonOrthoCorrector>
void Diffusion<NonOrthoCorrector>::apply_interior(const mesh::Face& face) {
    const auto& owner = _mesh.cell(face.owner());
    const auto& neighbor = _mesh.cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    const auto& [Sf, Ef, Tf] = _corrector.interior_triplet(owner, neighbor, face);

    // geometric diffusion coefficient
    auto g_diff = Ef.norm() / (d_CF_norm + EPSILON);

    auto owner_id = owner.id();
    auto neighbor_id = neighbor.id();

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    auto dval = g_diff * _kappa;
    insert(owner_id, owner_id, dval);
    insert(neighbor_id, neighbor_id, dval);

    // off-diagonal coefficients
    insert(owner_id, neighbor_id, -dval);
    insert(neighbor_id, owner_id, -dval);

    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukallad et al., 2015)
    auto grad_f = _corrector.grad_scheme().gradient_at_face(face);
    rhs(owner_id) += Tf.dot(grad_f) * _kappa;
    rhs(neighbor_id) += -Tf.dot(grad_f) * _kappa;
}

template <>
void inline Diffusion<nonortho::NilCorrector>::apply_interior(const mesh::Face& face) {
    const auto& owner = _mesh.cell(face.owner());
    const auto& neighbor = _mesh.cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    auto Sf = face.area_vector();

    // geometric diffusion coefficient
    auto g_diff = Sf.norm() / (d_CF_norm + EPSILON);

    auto owner_id = owner.id();
    auto neighbor_id = neighbor.id();

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    auto dval = g_diff * _kappa;
    insert(owner_id, owner_id, dval);
    insert(neighbor_id, neighbor_id, dval);

    // off-diagonal coefficients
    insert(owner_id, neighbor_id, -dval);
    insert(neighbor_id, owner_id, -dval);
}


template <typename NonOrthoCorrector>
void Diffusion<NonOrthoCorrector>::correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                                                    const mesh::Face& face,
                                                                    const Vector3d& T_f) {
    // we need to calculate the gradient of phi at the face
    // first let's calculate the vector joining the face center to the cell center
    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();
    auto e = d_CF / d_CF_norm;

    // now we need to calculate the gradient of phi at the face center
    auto boundary_patch_id = face.boundary_patch_id().value();
    const auto& face_boundary_patch = _mesh.boundary_patches()[boundary_patch_id];

    auto phi_wall = face_boundary_patch.get_scalar_bc(_phi.name());
    auto phi_c = _phi[cell.id()];

    auto grad_f = ((phi_wall - phi_c) / (d_CF_norm + EPSILON)) * e;
    rhs(cell.id()) += T_f.dot(grad_f) * _kappa;
}


template <typename NonOrthoCorrector>
void Diffusion<NonOrthoCorrector>::apply_boundary_fixed(const mesh::Cell& cell,
                                                        const mesh::Face& face) {
    // get the fixed phi variable associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - cell.center();
    auto d_Cf_norm = d_Cf.norm();

    const auto& [_, Ef, Tf] = _corrector.boundary_triplet(cell, face);
    auto g_diff = Ef.norm() / (d_Cf_norm + EPSILON);

    insert(cell_id, cell_id, g_diff * _kappa);
    rhs(cell_id) += g_diff * _kappa * phi_wall;

    correct_non_orhto_boundary_fixed(cell, face, Tf);
}

template <>
auto inline Diffusion<nonortho::NilCorrector>::requires_correction() const -> bool {
    return false;
}


} // namespace prism::diffusion