#include <algorithm>

#include "gradient.h"
#include "prism/constants.h"

namespace prism::gradient {

LeastSquares::LeastSquares(field::IScalar* field) : IGradient(field) {
    const auto& mesh = this->field()->mesh();
    _cell_gradients.reserve(mesh.nCells());

    std::transform(_cell_gradients.begin(),
                   _cell_gradients.end(),
                   std::back_inserter(_cell_gradients),
                   [](const auto& _) { return Vector3d::Zero(); }); // NOLINT

    setPseudoInvMatrices();
}

void LeastSquares::setPseudoInvMatrices() {
    // This function is based on section 9.3 'Least-Square Gradient'
    const auto& mesh = this->field()->mesh();

    // resize the pseudo-inverse matrices vector
    _pinv_matrices.resize(mesh.nCells());

    for (const auto& cell : mesh.cells()) {
        // A 3x3 matrix of the left hand side of equation (9.27)
        // for the k-th cell, we calculate this distance matrix D
        // and push the pseudo-inverse [(D * D^T)^{-1} * D^T] to _pinv_matrix vector
        Matrix3d d_matrix = Matrix3d::Zero();

        for (auto face_id : cell.facesIds()) {
            const auto& face = mesh.face(face_id);

            // This will hold the distance vector from neighbor cell center to k-th cell
            // center, or in case we have a boundary face, r_CF will be the distance vector
            // from boundary face center to the k-th cell center.
            // check equation (9.22)
            Vector3d r_CF = {.0, .0, .0};

            if (face.isInterior()) {
                // interior face
                const auto neighbor = mesh.otherSharingCell(cell, face);
                r_CF = neighbor.center() - cell.center();
            } else {
                // boundary face
                r_CF = face.center() - cell.center();
            }

            // weight factor defined in equation (9.28)
            const double wk = 1 / (r_CF.norm() + EPSILON);
            const double dx = r_CF.x(); // equation (9.24)
            const double dy = r_CF.y(); // equation (9.24)
            const double dz = r_CF.z(); // equation (9.24)

            // clang-format off
            // left hand side of equation (9.27) for the k-th cell, before summing and before
            // multiplying the weight factor wk
            Matrix3d d_matrix_k = Matrix3d::Zero();
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

auto LeastSquares::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    const auto& mesh = this->field()->mesh();

    // right hand side of equation (9.27)
    Vector3d b {0.0, 0.0, 0.0};

    for (auto face_id : cell.facesIds()) {
        const auto& face = mesh.face(face_id);

        double delta_phi = 0.0;
        auto phi_cell = this->field()->valueAtCell(cell);
        Vector3d r_CF = {.0, .0, .0};

        if (face.isInterior()) {
            // interior face
            const auto neighbor = mesh.otherSharingCell(cell, face);
            r_CF = neighbor.center() - cell.center();
            auto nei_phi = this->field()->valueAtCell(neighbor);
            delta_phi = nei_phi - phi_cell;

        } else {
            // boundary face
            auto bface_phi = this->field()->valueAtFace(face);
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
    Vector3d grad = _pinv_matrices[cell.id()] * b;
    _cell_gradients[cell.id()] = grad;

    return grad;
}

auto LeastSquares::gradAtCellStored(const mesh::Cell& cell) -> Vector3d {
    return _cell_gradients[cell.id()];
}
} // namespace prism::gradient