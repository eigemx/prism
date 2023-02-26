#include "fvscheme.h"

namespace prism {
void FVScheme::apply(const mesh::Cell& cell, const mesh::Face& face) const {
    if (!face.has_neighbor()) {
        apply_boundary(cell, face);
    } else {
        apply_interior(cell, face);
    }
}

} // namespace prism