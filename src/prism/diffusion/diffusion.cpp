#include "diffusion.h"

#include <cmath>

#include "../print.h"

namespace prism::diffusion {

Linear::Linear(double kappa, ScalarField& phi) : _kappa(kappa), _phi(phi), _mesh(phi.mesh()) {
    init_linear_system(_mesh);
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
    auto S_f = face.area() * face.normal();

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
    auto S_f = face.area() * face.normal();
    auto S_f_norm = face.area();

    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();

    // diffusion factor
    auto g_diff = (_kappa * S_f_norm) / (d_CF_norm + 1e-8);

    coeff_matrix().coeffRef(cell_id, cell_id) += g_diff;
    rhs_vector()[cell_id] += g_diff * phi_wall;
}
} // namespace prism::diffusion