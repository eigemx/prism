#include "equation.h"

namespace prism {
ConservedScalarSteadyEquation::ConservedScalarSteadyEquation(std::string scalar_name,
                                                             mesh::PMesh& mesh)
    : _scalar_name(std::move(scalar_name)), _mesh(mesh) {
    // reserve space for the coefficient matrix with the size of mesh cells
    _coeff_matrix.reserve(mesh.cells().size());

    // reserve space for the right hand side vector with the size of mesh cells
    _b.resize(mesh.cells().size());
}
} // namespace prism
