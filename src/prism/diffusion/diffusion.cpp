#include "diffusion.h"

#include <cmath>

#include "../print.h"

namespace prism::diffusion {

void Linear::apply_interior(const mesh::Cell& cell,
                            const mesh::Face& face,
                            const mesh::PMesh& mesh,
                            SparseMatrix& coeffs_matrix,
                            VectorXd& rhs_vec) const {
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

    const auto& adj_cell = mesh.cells()[adj_cell_id];

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
    coeffs_matrix.coeffRef(cell_id, cell_id) += g_diff;
    coeffs_matrix.coeffRef(cell_id, adj_cell_id) += -g_diff;
}

void Linear::apply_boundary(const mesh::Cell& cell,
                            const mesh::Face& face,
                            const mesh::PMesh& mesh,
                            SparseMatrix& coeffs_matrix,
                            VectorXd& rhs_vec) const {
    auto cell_id = cell.id();
    const auto& boundary_patch = mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        case mesh::BoundaryPatchType::Wall: {
            apply_wall_boundary(cell, face, mesh, coeffs_matrix, rhs_vec);
            return;
        }

        default:
            throw std::runtime_error("Unsupported boundary type");
    }
}

void Linear::apply_wall_boundary(const mesh::Cell& cell,
                                 const mesh::Face& face,
                                 const mesh::PMesh& mesh,
                                 SparseMatrix& coeffs_matrix,
                                 VectorXd& rhs_vec) const {
    const auto& boundary_patch = mesh.face_boundary_patch(face);

    const auto& wall_bdata = std::get<mesh::WallBoundaryData>(boundary_patch.data());
    auto phi_wall = wall_bdata.temperature.value();

    auto cell_id = cell.id();

    // face area vector
    auto S_f = face.area() * face.normal();
    auto S_f_norm = face.area();

    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();

    // diffusion factor
    auto g_diff = (_kappa * S_f_norm) / (d_CF_norm + 1e-8);

    coeffs_matrix.coeffRef(cell_id, cell_id) += g_diff;
    rhs_vec[cell_id] += g_diff * phi_wall;
}
} // namespace prism::diffusion