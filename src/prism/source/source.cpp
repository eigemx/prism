#include "source.h"

namespace prism::source {
ConstantScalar::ConstantScalar(ScalarField& phi) : _phi(phi), FVScheme(phi.mesh().n_cells()) {
    _volume_field.resize(phi.mesh().n_cells());

    for (const auto& cell : phi.mesh().cells()) {
        _volume_field[cell.id()] = cell.volume();
    }
}

void ConstantScalar::apply() {
    rhs() = _phi.data().array() * _volume_field.array();
}

// TODO: replace calls to apply interior by simply making coeff_matrix() an I matrix,
// or -I matrix in case of SourceSign::Positive, during construction of ImplicitPhi

//template <>
//void ImplicitPhi<SourceSign::Positive>::apply_interior(const mesh::Cell& cell,
//                                                       const mesh::Face& face) //NOLINT
//{
//    auto cell_id = cell.id();
//    matrix(cell_id, cell_id) = -_coeff;
//}
//
//template <>
//void ImplicitPhi<SourceSign::Negative>::apply_interior(const mesh::Cell& cell,
//                                                       const mesh::Face& face) // NOLINT
//{
//    auto cell_id = cell.id();
//    matrix(cell_id, cell_id) = _coeff;
//}

} // namespace prism::source