#include "equation.h"

#include "print.h"

namespace prism {
SteadyConservedScalar::SteadyConservedScalar(std::string scalar_name, mesh::PMesh& mesh)
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
                auto altered_coeffs = scheme->apply(cell, face, _mesh);

                // update cell coefficient
                coeff_matrix().coeffRef(cell.id(), cell.id()) += altered_coeffs.central;

                // update neighbour coefficients
                if (altered_coeffs.neighbor.has_value()) {
                    auto nei_cell_id = altered_coeffs.neighbor.value().first;
                    auto coeff = altered_coeffs.neighbor.value().second;
                    coeff_matrix().coeffRef(cell.id(), nei_cell_id) += coeff;
                }

                // update right hand side vector
                lhs_vector()(cell.id()) += altered_coeffs.b;
            }
        }
    }
}
} // namespace prism
