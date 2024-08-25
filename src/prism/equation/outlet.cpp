#include "boundary.h"

namespace prism::eqn::boundary {
void Outlet<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {}

} // namespace prism::eqn::boundary