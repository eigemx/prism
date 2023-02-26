#pragma once

#include <memory>
#include <vector>

#include "equation.h"
#include "fvscheme.h"

namespace prism {
class CellVisitor {
  public:
    CellVisitor(const mesh::PMesh& mesh,
                const std::vector<std::shared_ptr<FVScheme>>& schemes,
                const LinearSystem& system)
        : _mesh(mesh), _schemes(schemes), _system(system) {}

    void visit_cells();

  private:
    const std::vector<std::shared_ptr<FVScheme>>& _schemes;
    const mesh::PMesh& _mesh;
    const LinearSystem& _system;
};
} // namespace prism
