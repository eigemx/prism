#include "source.h"

namespace prism::source {
// TODO: replace calls to apply interior by simply making coeff_matrix() an I matrix,
// or -I matrix in case of SourceSign::Positive, during construction of ImplicitPhi

template <>
void ImplicitPhi<SourceSign::Positive>::apply_interior(const mesh::Cell& cell,
                                                       const mesh::Face& face) //NOLINT
{
    auto cell_id = cell.id();
    coeff_matrix().coeffRef(cell_id, cell_id) = -_coeff;
}

template <>
void ImplicitPhi<SourceSign::Negative>::apply_interior(const mesh::Cell& cell,
                                                       const mesh::Face& face) // NOLINT
{
    auto cell_id = cell.id();
    coeff_matrix().coeffRef(cell_id, cell_id) = _coeff;
}

} // namespace prism::source