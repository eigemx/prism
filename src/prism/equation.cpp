#include "equation.h"

#include "print.h"

namespace prism {

Equation::Equation(ScalarField& phi, std::vector<FVScheme*> schemes)
    : _schemes(std::move(schemes)), _phi(phi) {
    if (_schemes.empty()) {
        throw std::runtime_error(
            format("Equation constructor was called with an empty schemes vector. No matrix "
                   "initialization was performed."));
    }

    const auto& mesh = _phi.mesh();
    auto n_cells = mesh.cells().size();

    _unified_coeff_matrix = SparseMatrix(n_cells, n_cells);
    _unified_rhs_vector = VectorXd::Zero(n_cells);
}

void Equation::update_coeffs() {
    if (_schemes.empty()) {
        throw std::runtime_error(
            "Equation::update_coeffs(): Attempting to update matrix coefficients with no FV "
            "schemes");
    }

    const auto& mesh = _schemes[0]->mesh();

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

    // zero out the matrix and vector for the linear system
    reset();

    for (auto* scheme : _schemes) {
        _unified_coeff_matrix += scheme->coeff_matrix();
        _unified_rhs_vector += scheme->rhs_vector();

        // prepare the scheme for the next iteration
        scheme->finalize();
    }
}

void Equation::reset() {
    _unified_coeff_matrix.setZero();
    _unified_rhs_vector.setZero();
}

} // namespace prism
