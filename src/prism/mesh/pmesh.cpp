#include "pmesh.h"

namespace prism::mesh {
// define all functions in PMesh class in pmesh.h
PMesh::PMesh(std::vector<Cell>&& cells, std::vector<Face>&& faces)
    : _cells(std::move(cells)), _faces(std::move(faces)) {}

PMesh::PMesh(std::vector<Cell>&& cells, std::vector<Face>&& faces,
             BoundaryConditions&& boundary_patches)
    : _cells(std::move(cells)),
      _faces(std::move(faces)),
      _boundary_conditions(std::move(boundary_patches)) {}

void PMesh::set_boundary_conditions(BoundaryConditions&& boundary_patches) {
    _boundary_conditions = std::move(boundary_patches);
}

void PMesh::set_boundary_conditions(const BoundaryConditions& boundary_patches) {
    _boundary_conditions = boundary_patches;
}

} // namespace prism::mesh