#include "visitor.h"

#include "print.h"

namespace prism {
void CellVisitor::visit_cells() {
    for (const auto& cell : _mesh.cells()) {
        for (const auto& face_id : cell.faces_ids()) {
            const auto& face = _mesh.faces()[face_id];

            for (auto& scheme : _schemes) {
                if (scheme->run_once() && scheme->is_first_run_complete()) {
                    continue;
                }
                scheme->apply(cell, face);
            }
        }
    }

    // first run is now complete, we need to mark schemes that require just one run as such
    for (auto& scheme : _schemes) {
        if (scheme->run_once()) {
            scheme->set_first_run_complete();
        }
    }
}
} // namespace prism