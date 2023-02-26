#include "visitor.h"

#include "print.h"

namespace prism {
void CellVisitor::visit_cells() {
    fmt::print("Visiting cells...\n");
    for (const auto& cell : _mesh.cells()) {
        for (const auto& face_id : cell.faces_ids()) {
            const auto& face = _mesh.faces()[face_id];

            for (const auto& scheme : _schemes) {
                scheme->apply(cell, face);
            }
        }
    }
}
} // namespace prism