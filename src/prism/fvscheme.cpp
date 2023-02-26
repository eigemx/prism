#include "fvscheme.h"

#include "print.h"

namespace prism {
void FVScheme::apply(const mesh::Cell& cell, const mesh::Face& face) const {
    fmt::print("Applying scheme to cell {} and face {}\n", cell.id(), face.id());
    /*if (!face.has_neighbor()) {
        apply_boundary(cell, face);
    } else {
        apply_interior(cell, face);
    }*/
}

} // namespace prism