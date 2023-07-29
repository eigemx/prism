#include "equation.h"

#include "print.h"

namespace prism {

void Equation::update_coeffs() {
    const auto& mesh = _phi.mesh();

    // iterate over all cells
    for (const auto& cell : mesh.cells()) {
        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.face(face_id);

            // iterate over all equation's finite volume schemes
            for (auto& scheme : _schemes) {
                // apply the scheme to the cell and face
                scheme->apply(cell, face);
            }
        }
    }

    for (auto& scheme : _schemes) {
        // prepare the scheme for the next iteration
        scheme->finalize();
    }

    // update the universal coefficient matrix and RHS vector
    for (auto& scheme : _schemes) {
        matrix() += scheme->matrix();
        rhs() += scheme->rhs();
    }
}

void Equation::zero_out_coeffs() {
    for (auto& scheme : _schemes) {
        if (scheme->requires_correction()) {
            scheme->matrix().setZero();
            scheme->rhs().setZero();
        }
    }

    // zero out the universal coefficient matrix and RHS vector
    matrix().setZero();
    rhs().setZero();
}

} // namespace prism
