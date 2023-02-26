#include "diffusion.h"

namespace prism {

void LinearDiffusionScheme::apply_interior(const mesh::Cell& cell, const mesh::Face& face) const {
    // apply orthogonal-like diffusion
    auto cell_id = cell.id();
    auto is_owned = face.owner() == cell_id;
}

void LinearDiffusionScheme::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const {
    // apply orthogonal-like diffusion
    auto cell_id = cell.id();
    const auto& boundary_patch = _mesh.face_boundary_patch(face); // TODO: This is awfully slow

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty:
            break;

        case mesh::BoundaryPatchType::Wall:
            break;
    }
}
} // namespace prism