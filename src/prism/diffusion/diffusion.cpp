#include "diffusion.h"

#include <cmath>

#include "../print.h"

namespace prism::diffusion {

Linear::Linear(double kappa, ScalarField& phi)
    : _kappa(kappa),
      _phi(phi),
      _mesh(phi.mesh()),
      _gradient_scheme(std::make_unique<gradient::GreenGauss>()) {
    init_linear_system(_mesh);
}

template <typename GradientScheme>
Linear::Linear(double kappa, ScalarField& phi, const GradientScheme& gradient_scheme)
    : _kappa(kappa),
      _phi(phi),
      _mesh(phi.mesh()),
      _gradient_scheme(std::make_unique<GradientScheme>(gradient_scheme)) {
    init_linear_system(_mesh);
}

void Linear::finalize() {
    // zero out the linear system to prepare for the next iteration
    coeff_matrix().setZero();
    rhs_vector().setZero();
}


void Linear::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();

    // get adjacent cell sharing `face` with `cell`
    std::size_t adj_cell_id {};

    if (face.owner() == cell_id) {
        // `face` is owned by `cell`, get the neighbor cell id
        adj_cell_id = face.neighbor().value();
    } else {
        // `face` is a neighbor to `cell`, get its owner id
        adj_cell_id = face.owner();
    }

    const auto& adj_cell = _mesh.cells()[adj_cell_id];

    // face area vector
    const auto& S_f = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = adj_cell.center() - cell.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    // orthogonal-like normal vector E_f using over-relaxed approach
    auto E_f = ((S_f.dot(S_f) / (e.dot(S_f) + 1e-8))) * e;

    // diffusion factor
    auto g_diff = (E_f.norm() * _kappa) / (d_CF_norm + 1e-8);

    // g_diff * (Φ_c - Φ_f)
    coeff_matrix().coeffRef(cell_id, cell_id) += g_diff;
    coeff_matrix().coeffRef(cell_id, adj_cell_id) += -g_diff;

    correct_non_orhto_interior(cell, adj_cell, face, S_f - E_f);
}


void Linear::correct_non_orhto_interior(const mesh::Cell& cell,
                                        const mesh::Cell& nei_cell,
                                        const mesh::Face& face,
                                        const Vector3d& T_f) {
    auto grad_c = _gradient_scheme->gradient(cell, _phi);
    auto grad_n = _gradient_scheme->gradient(nei_cell, _phi);
    auto gc = mesh::PMesh::cells_weighting_factor(cell, nei_cell, face);
    auto grad_f = (gc * grad_c) + ((1 - gc) * grad_n);

    rhs_vector()[cell.id()] += T_f.dot(grad_f) * _kappa;
}


void Linear::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();
    const auto& boundary_patch = _mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        case mesh::BoundaryPatchType::Fixed: {
            apply_boundary_fixed(cell, face);
            return;
        }

        default:
            throw std::runtime_error("Unsupported boundary type");
    }
}


void Linear::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_attribute(_phi.name());
    auto cell_id = cell.id();

    // face area vector
    const auto& S_f = face.area_vector();
    auto S_f_norm = face.area();

    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();
    auto e = d_CF / d_CF_norm;

    auto E_f = ((S_f.dot(S_f) / (e.dot(S_f) + 1e-8))) * e;

    // diffusion factor
    auto g_diff = (E_f.norm() * _kappa) / (d_CF_norm + 1e-8);

    coeff_matrix().coeffRef(cell_id, cell_id) += g_diff;
    rhs_vector()[cell_id] += g_diff * phi_wall;

    correct_non_orhto_boundary_fixed(cell, face, S_f - E_f);
}


void Linear::correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                              const mesh::Face& face,
                                              const Vector3d& T_f) {
    // we need to calculate the gradient of phi at the face
    // first let's calculate the vector joining the face center to the cell center
    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();
    auto e = d_CF / d_CF_norm;

    // now we need to calculate the gradient of phi at the face center
    const auto& face_boundary_patch = _mesh.boundary_patches()[face.boundary_patch_id().value()];
    auto phi_wall = face_boundary_patch.get_scalar_attribute(_phi.name());
    auto phi_c = _phi[cell.id()];

    auto grad_f = ((phi_wall - phi_c) / d_CF_norm) * e;

    rhs_vector()[cell.id()] += T_f.dot(grad_f) * _kappa;
}
} // namespace prism::diffusion