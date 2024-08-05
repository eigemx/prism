#include <spdlog/spdlog.h>

#include "gradient.h"
#include "prism/constants.h"
#include "prism/types.h"

namespace prism::gradient {
LeastSquares::LeastSquares(const field::Scalar& field) : _field(field), IGradient(field) {
    set_pseudo_inv_matrices();
}

void LeastSquares::set_pseudo_inv_matrices() {
    // This function is based on section 9.3 'Least-Square Gradient'
    const auto& mesh = _field.mesh();

    // resize the pseudo-inverse matrices vector
    _pinv_matrices.resize(mesh.n_cells());

    for (const auto& cell : mesh.cells()) {
        // A 3x3 matrix of the left hand side of equation (9.27)
        // for the k-th cell, we calculate this distance matrix D
        // and push the pseudo-inverse [(D * D^T)^{-1} * D^T] to _pinv_matrix vector
        Matrix3d d_matrix = Matrix3d::Zero();

        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.face(face_id);

            // This will hold the distance vector from neighbor cell center to k-th cell
            // center, or in case we have a boundary face, r_CF will be the distance vector
            // from boundary face center to the k-th cell center.
            // check equation (9.22)
            Vector3d r_CF = {.0, .0, .0};

            // difference of field values between the k-th cell and its i-th neigbor
            // or its boundary face field value
            double delta_phi = 0.0;

            Matrix3d d_matrix_k = Matrix3d::Zero();

            if (face.is_interior()) {
                // interior face
                const auto neighbor = mesh.other_sharing_cell(cell, face);
                r_CF = neighbor.center() - cell.center();
                delta_phi = _field.valueAtCell(neighbor) - _field.valueAtCell(cell);
            } else {
                // boundary face
                r_CF = face.center() - cell.center();
                delta_phi = _field.valueAtFace(face) - _field.valueAtCell(cell);
            }

            // weight factor defined in equation (9.28)
            const double wk = 1 / (r_CF.norm() + EPSILON);
            const double dx = r_CF.x(); // equation (9.24)
            const double dy = r_CF.y(); // equation (9.24)
            const double dz = r_CF.z(); // equation (9.24)

            // clang-format off
            // left hand side of equation (9.27) for the k-th cell, before summing and before
            // multiplying the weight factor wk
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
    const auto& mesh = _field.mesh();

    // right hand side of equation (9.27)
    Vector3d b {0.0, 0.0, 0.0};

    for (auto face_id : cell.faces_ids()) {
        const auto& face = mesh.face(face_id);

        double delta_phi = 0.0;
        auto phi_cell = _field.valueAtCell(cell);
        Vector3d r_CF = {.0, .0, .0};

        if (face.is_interior()) {
            // interior face
            const auto neighbor = mesh.other_sharing_cell(cell, face);
            r_CF = neighbor.center() - cell.center();
            auto nei_phi = _field.valueAtCell(neighbor);
            delta_phi = nei_phi - phi_cell;

        } else {
            // boundary face
            auto bface_phi = _field.valueAtFace(face);
            r_CF = face.center() - cell.center();
            delta_phi = bface_phi - phi_cell;
        }
        const double wk = 1 / (r_CF.norm() + EPSILON);
        const double dx = r_CF.x();
        const double dy = r_CF.y();
        const double dz = r_CF.z();

        // update the right hand side of equation (9.27)
        b += Vector3d {dx * delta_phi, dy * delta_phi, dz * delta_phi} * wk;
    }
    return _pinv_matrices[cell.id()] * b;
}

} // namespace prism::gradient