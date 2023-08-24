#include "diffusion.h"

#include <cmath>

#include "../print.h"
#include "../types.h"

namespace prism::diffusion {

template <>
void Diffusion<NonOrthoCorrection::OverRelaxed>::apply_interior(const mesh::Face& face) {
    const auto& owner = _mesh.cell(face.owner());
    const auto& neighbor = _mesh.cell(face.neighbor().value());
    const auto& S_f = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    // orthogonal-like normal vector E_f using over-relaxed approach
    auto E_f = ((S_f.dot(S_f) / (e.dot(S_f) + EPSILON))) * e;

    // geometric diffusion coefficient
    auto g_diff = E_f.norm() / (d_CF_norm + EPSILON);

    auto owner_id = owner.id();
    auto neighbor_id = neighbor.id();

    // kappa * g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    matrix(owner_id, owner_id) += g_diff * _kappa;
    matrix(neighbor_id, neighbor_id) += g_diff * _kappa;

    // off-diagonal coefficients
    matrix(owner_id, neighbor_id) += -g_diff * _kappa;
    matrix(neighbor_id, owner_id) += -g_diff * _kappa;

    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukallad et al., 2015)
    auto grad_f = _gradient_scheme->gradient_at_face(face);
    auto T_f = S_f - E_f;
    rhs(owner_id) += T_f.dot(grad_f) * _kappa;
    rhs(neighbor_id) += -T_f.dot(grad_f) * _kappa;
}

template <>
void Diffusion<NonOrthoCorrection::OverRelaxed>::correct_non_orhto_boundary_fixed(
    const mesh::Cell& cell,
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

template <>
void Diffusion<NonOrthoCorrection::OverRelaxed>::apply_boundary_fixed(const mesh::Cell& cell,
                                                                      const mesh::Face& face) {
    // get the fixed phi variable associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    const auto& S_f = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - cell.center();
    auto d_Cf_norm = d_Cf.norm();
    auto e = d_Cf / d_Cf_norm;
    auto E_f = ((S_f.dot(S_f) / e.dot(S_f))) * e;

    auto g_diff = E_f.norm() / (d_Cf_norm + EPSILON);

    matrix(cell_id, cell_id) += g_diff * _kappa;
    rhs(cell_id) += g_diff * _kappa * phi_wall;

    correct_non_orhto_boundary_fixed(cell, face, S_f - E_f);
}

template <>
auto Diffusion<NonOrthoCorrection::OverRelaxed>::requires_correction() const -> bool {
    return true;
}

template <>
void Diffusion<NonOrthoCorrection::None>::apply_interior(const mesh::Face& face) {
    const auto& owner = _mesh.cell(face.owner());
    const auto& neighbor = _mesh.cell(face.neighbor().value());

    auto owner_id = owner.id();
    auto neighbor_id = neighbor.id();

    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // geometric diffusion coefficient
    auto g_diff = face.area() / (d_CF_norm + EPSILON);

    // kappa * g_diff * (Φ_C - Φ_N)
    matrix(owner_id, owner_id) += g_diff * _kappa;
    matrix(neighbor_id, neighbor_id) += g_diff * _kappa;

    matrix(owner_id, neighbor_id) += -g_diff * _kappa;
    matrix(neighbor_id, owner_id) += -g_diff * _kappa;
}

template <>
void Diffusion<NonOrthoCorrection::None>::apply_boundary_fixed(const mesh::Cell& cell,
                                                               const mesh::Face& face) {
    // get the fixed phi variable associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    const auto& S_f = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - cell.center();
    auto d_Cf_norm = d_Cf.norm();

    auto g_diff = face.area() / (d_Cf_norm + EPSILON);

    matrix(cell_id, cell_id) += g_diff * _kappa;
    rhs(cell_id) += g_diff * _kappa * phi_wall;
}

template <>
auto Diffusion<NonOrthoCorrection::None>::requires_correction() const -> bool {
    return false;
}

} // namespace prism::diffusion