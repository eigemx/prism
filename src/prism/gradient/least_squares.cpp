#include "gradient.h"
#include "prism/constants.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"

namespace prism::gradient {
LeastSquares::LeastSquares(const ScalarField& field) : _field(field), AbstractGradient(field) {
    set_lsq_matrices();
}

void LeastSquares::set_lsq_matrices() {
    const auto& mesh = _field.mesh();
    _pinv_matrices.resize(mesh.n_cells());

    for (const auto& cell : mesh.cells()) {
        Matrix3d d_matrix = Matrix3d::Zero();

        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.face(face_id);

            Vector3d r_CF = {.0, .0, .0};
            double delta_phi = 0.0;
            Matrix3d d_matrix_k = Matrix3d::Zero();

            if (face.is_interior()) {
                // interior face
                const auto neighbor = mesh.other_sharing_cell(cell, face);
                r_CF = neighbor.center() - cell.center();
                delta_phi = _field.value_at_cell(neighbor) - _field.value_at_cell(cell);
            } else {
                // boundary face
                r_CF = face.center() - cell.center();
                delta_phi = _field.value_at_face(face) - _field.value_at_cell(cell);
            }

            const double wk = 1 / (r_CF.norm() + EPSILON);
            const double dx = r_CF.x();
            const double dy = r_CF.y();
            const double dz = r_CF.z();

            // clang-format off
            d_matrix_k << (dx * dx), (dx * dy), (dx * dz), 
                          (dy * dx), (dy * dy), (dy * dz),
                          (dz * dx), (dz * dy), (dz * dz);
            // clang-format on

            d_matrix += d_matrix_k * wk;
        }

        const auto& d_matrix_t = d_matrix.transpose();
        _pinv_matrices[cell.id()] = (d_matrix_t * d_matrix).inverse() * d_matrix_t;
    }
}

auto LeastSquares::gradient_at_cell(const mesh::Cell& cell) -> Vector3d {
    const auto n_faces = cell.faces_ids().size();
    const auto& mesh = _field.mesh();
    const auto& inverse_matrix = _pinv_matrices[cell.id()];

    Vector3d b {0.0, 0.0, 0.0};

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.face(face_id);

        double delta_phi = 0.0;
        auto phi_cell = _field.value_at_cell(cell);
        Vector3d r_CF = {.0, .0, .0};

        if (face.is_interior()) {
            // interior face
            const auto neighbor = mesh.other_sharing_cell(cell, face);
            r_CF = neighbor.center() - cell.center();
            auto nei_phi = _field.value_at_cell(neighbor);
            delta_phi = nei_phi - phi_cell;

        } else {
            // boundary face
            auto bface_phi = _field.value_at_face(face);
            r_CF = face.center() - cell.center();
            delta_phi = bface_phi - phi_cell;
        }
        const double wk = 1 / (r_CF.norm() + EPSILON);
        const double dx = r_CF.x();
        const double dy = r_CF.y();
        const double dz = r_CF.z();
        b += Vector3d {dx * delta_phi, dy * delta_phi, dz * delta_phi} * wk;
    }
    return inverse_matrix * b;
}

} // namespace prism::gradient