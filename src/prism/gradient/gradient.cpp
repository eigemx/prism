#include "gradient.h"

namespace prism::gradient {
void GreenGauss::apply_interior(const mesh::Cell& cell, const mesh::Face& face) const {}
void GreenGauss::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const {}

} // namespace prism::gradient