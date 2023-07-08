#include "source.h"

namespace prism::source {

// NOLINTNEXTLINE (face is unused, but required by interface)
void ConstantScalar::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto q = _phi[cell.id()];
    rhs_vector()(cell.id()) = q * cell.volume();
}
} // namespace prism::source