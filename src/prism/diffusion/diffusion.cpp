#include "diffusion.h"

#include <cmath>

#include "../print.h"

namespace prism::diffusion {

void Linear::apply_interior(const mesh::Cell& cell, const mesh::Face& face) const {
    auto cell_id = cell.id();

    // does `cell` own `face` or it is a neighbor?
    auto is_owned = (face.owner() == cell_id);
    auto own_factor = is_owned ? 0.0 : 1.0;

    // get adjacent cell shaing `face` with `cell`
    auto adj_cell_id = is_owned ? face.neighbor().value() : face.owner();
    auto adj_cell = _mesh.cells()[adj_cell_id];

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
    auto g_diff = (-E_f.norm() * _kappa) / (d_CF_norm + 1e-8);
    g_diff *= std::pow(-1., own_factor);

    // update coefficients matrix
    _coeffs.coeffRef(cell_id, cell_id) += -g_diff;
    _coeffs.coeffRef(cell_id, adj_cell_id) += g_diff;
}

void Linear::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const {
    auto cell_id = cell.id();
    const auto& boundary_patch = _mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty:
            break;

        case mesh::BoundaryPatchType::Wall:
            auto T_wall =
                std::get<mesh::WallBoundaryData>(boundary_patch.data()).temperature.value();
            apply_wall_boundary(cell, face, T_wall);
    }
}

void Linear::apply_wall_boundary(const mesh::Cell& cell,
                                 const mesh::Face& face,
                                 double phi_wall) const {
    auto cell_id = cell.id();

    // face area vector
    auto S_f = face.area() * face.normal();
    auto S_f_norm = face.area();

    auto d_CF = face.center() - cell.center();
    auto d_CF_norm = d_CF.norm();

    // diffusion factor
    auto g_diff = -_kappa * S_f_norm / (d_CF_norm + 1e-8);

    // update coefficients matrix
    _coeffs.coeffRef(cell_id, cell_id) += -g_diff;
    _b[cell_id] += -g_diff * phi_wall;
}
} // namespace prism::diffusion