#include "../mesh/utilities.h"
#include "gradient.h"

// TODO: Validate the results of using LeastSquares as an explicit gradient calculator,
// also check if skewness correction is required here or not, as face_grad is being calculated
// by weighting the gradients of the two adjacent cells.
namespace prism::gradient {
LeastSquares::LeastSquares(const ScalarField& field) : _field(field), AbstractGradient(field) {
    set_lsq_matrices();
}

void LeastSquares::set_lsq_matrices() {
    const auto& mesh = _field.mesh();
    _pinv_matrices.resize(mesh.n_cells());

    for (const auto& cell : mesh.cells()) {
        auto c_id = cell.id();
        auto n_faces = cell.faces_ids().size();

        MatrixX3d d_matrix = MatrixX3d::Zero(n_faces, 3);

        for (std::size_t i = 0; i < n_faces; ++i) {
            auto f_id = cell.faces_ids()[i];
            const auto& face = mesh.face(f_id);

            Vector3d d_CF {.0, .0, .0};

            if (face.has_neighbor()) {
                // interior face
                const auto neighbor = mesh.other_sharing_cell(cell, face);
                d_CF = neighbor.center() - cell.center();
            } else {
                // boundary face
                d_CF = face.center() - cell.center();
            }

            d_matrix.row(i) = d_CF;
        }

        const auto& d_matrix_t = d_matrix.transpose();
        _pinv_matrices[c_id] = (d_matrix_t * d_matrix).inverse() * d_matrix_t;
    }
}

auto LeastSquares::gradient_at_cell(const mesh::Cell& cell) -> Vector3d {
    auto c_id = cell.id();
    const auto n_faces = cell.faces_ids().size();

    const auto& mesh = _field.mesh();

    const auto& inverse_matrix = _pinv_matrices[c_id];

    VectorXd phi_diff = VectorXd::Zero(n_faces);

    Vector3d grad {.0, .0, .0};

    for (std::size_t i = 0; i < n_faces; ++i) {
        auto f_id = cell.faces_ids()[i];
        const auto& face = mesh.face(f_id);
        auto phi = _field[cell.id()];

        if (face.has_neighbor()) {
            // interior face
            const auto neighbor = mesh.other_sharing_cell(cell, face);
            auto nei_phi = _field[neighbor.id()];
            phi_diff[i] = nei_phi - phi;

        } else {
            // boundary face
            phi_diff[i] = boundary_face_phi(face) - phi;
        }
    }

    grad = inverse_matrix * phi_diff;

    return grad;
}

auto LeastSquares::boundary_face_phi(const mesh::Face& f) -> double {
    const auto& mesh = _field.mesh();
    const auto& patch = mesh.boundary_patch(f);
    const auto& bc = patch.get_bc(_field.name());

    switch (bc.bc_type()) {
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet:
            return patch.get_scalar_bc(_field.name());

        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Empty:
        case mesh::BoundaryConditionType::Outlet:
            return _field[f.owner()];

        default:
            throw std::runtime_error(
                fmt::format("gradient::LeastSquared::boundary_face_phi(): "
                            "Non-implemented boundary type for boundary patch: '{}'",
                            patch.name()));
    }
}

} // namespace prism::gradient