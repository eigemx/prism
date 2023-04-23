#include "equation.h"

#include "print.h"

namespace prism {

Equation::Equation(ScalarField& phi, std::vector<FVScheme*> schemes)
    : _schemes(std::move(schemes)), _phi(phi) {
    if (_schemes.empty()) {
        throw std::runtime_error(
            format("Equation constructor was called with an empty FV schemes vector. No matrix "
                   "initialization was performed."));
    }

    const auto& mesh = _phi.mesh();
    auto n_cells = mesh.cells().size();

    _coeff_matrix = SparseMatrix(n_cells, n_cells);
    _rhs_vector = VectorXd::Zero(n_cells);

    for (auto* scheme : _schemes) {
        scheme->connect_linear_system(_coeff_matrix, _rhs_vector);
    }
}

void Equation::update_coeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all cells
    for (const auto& cell : mesh.cells()) {
        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.faces()[face_id];

            // iterate over all equation finite volume schemes
            for (auto& scheme : _schemes) {
                // apply the scheme to the cell and face
                scheme->apply(cell, face);
            }
        }
    }

    for (auto* scheme : _schemes) {
        // prepare the scheme for the next iteration
        scheme->finalize();
    }
}

void Equation::reset_coeffs() {
    _coeff_matrix.setZero();
    _rhs_vector.setZero();
}

} // namespace prism
