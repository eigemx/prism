#include "diffusion.h"

namespace prism {

void LinearDiffusionScheme::apply_interior(const mesh::Cell& cell, const mesh::Face& face) const {
    auto cell_id = cell.id();
    auto adj_cell_id = face.neighbor();

    // does `cell` own `face` or it is a neighbor?
    auto is_owned = (face.owner() == cell_id);
}

void LinearDiffusionScheme::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const {
    auto cell_id = cell.id();
    const auto& boundary_patch = _mesh.face_boundary_patch(face);

    switch (boundary_patch.type()) {
        case mesh::BoundaryPatchType::Empty:
            break;

        case mesh::BoundaryPatchType::Wall:
            break;

        default:
            throw std::runtime_error("Non-implemented boundary patch type");
    }
}
} // namespace prism