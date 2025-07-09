
#include "ifield.h"

#include <stdexcept>

#include "prism/mesh/pmesh.h"

namespace prism::field {

IScalar::IScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh)
    : IField<double>(std::move(name), mesh) {}

IVector::IVector(std::string name, const SharedPtr<mesh::PMesh>& mesh)
    : IField<Vector3d>(std::move(name), mesh) {}
} // namespace prism::field

namespace prism::field::detail {

void checkFieldName(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

void checkMesh(const SharedPtr<mesh::PMesh>& mesh) {
    if (mesh->cells().empty() || mesh->faces().empty() || mesh->boundaryPatches().empty()) {
        throw std::runtime_error(
            "Cannot create a field over either an empty mesh or a mesh that have no boundary "
            "patches.");
    }
}
} // namespace prism::field::detail
