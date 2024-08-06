
#include "ifield.h"

#include <stdexcept>

#include "prism/mesh/pmesh.h"


namespace prism::field::detail {

void checkFieldName(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

void checkMesh(const mesh::PMesh& mesh) {
    if (mesh.cells().empty() || mesh.faces().empty() || mesh.boundaryPatches().empty()) {
        throw std::runtime_error("Cannot create a field over an empty mesh.");
    }
}

} // namespace prism::field::detail