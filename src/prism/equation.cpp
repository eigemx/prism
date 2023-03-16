#include "equation.h"

#include "print.h"

namespace prism {
SteadyConservedScalar::SteadyConservedScalar(std::string scalar_name, const mesh::PMesh& mesh)
    : _scalar_name(std::move(scalar_name)), _mesh(mesh) {
    // reserve space for the coefficient matrix with the size of mesh cells
    auto n_cells = mesh.cells().size();
    coeff_matrix() = SparseMatrix(n_cells, n_cells);

    // reserve space for the right hand side vector with the size of mesh cells
    lhs_vector().resize(n_cells);
    lhs_vector().setZero();
}

void SteadyConservedScalar::update_coeffs() {
    // iterate over all cells
    for (const auto& cell : _mesh.cells()) {
        for (auto face_id : cell.faces_ids()) {
            const auto& face = _mesh.faces()[face_id];

            // iterate over all schemes
            for (const auto& scheme : schemes()) {
                // apply the scheme to the cell and face
                scheme->apply(cell, face, _mesh, coeff_matrix(), lhs_vector());
            }
        }
    }
}
} // namespace prism
