#include "equation.h"

#include <cassert>

#include "print.h"

namespace prism {

Equation::Equation(ScalarField& phi, std::vector<FVScheme*> schemes)
    : _schemes(std::move(schemes)), _phi(phi) {
    if (_schemes.empty()) {
        throw std::runtime_error(
            fmt::format("Equation constructor was called with an empty FV schemes vector. "
                        "At least one scheme is required."));
    }

    const auto& mesh = _phi.mesh();
    auto n_cells = mesh.cells().size();

    assert(n_cells > 0 && "Equation constructor was called with an empty mesh.");

    _coeff_matrix = SparseMatrix(n_cells, n_cells);
    _rhs_vector = VectorXd::Zero(n_cells);
}

void Equation::update_coeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all cells
    for (const auto& cell : mesh.cells()) {
        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.face(face_id);

            // iterate over all equation finite volume schemes
            for (auto* scheme : _schemes) {
                // apply the scheme to the cell and face
                scheme->apply(cell, face);
            }
        }
    }

    for (auto* scheme : _schemes) {
        // prepare the scheme for the next iteration
        scheme->finalize();
    }

    // update the universal coefficient matrix and RHS vector
    for (auto* scheme : _schemes) {
        _coeff_matrix += scheme->coeff_matrix();
        _rhs_vector += scheme->rhs_vector();
    }
}

void Equation::zero_out_coeffs() {
    for (auto* scheme : _schemes) {
        if (scheme->requires_correction()) {
            scheme->coeff_matrix().setZero();
            scheme->rhs_vector().setZero();
        }
    }

    _coeff_matrix.setZero();
    _rhs_vector.setZero();
}

} // namespace prism
