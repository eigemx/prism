#pragma once

#include <memory>
#include <vector>

#include "equation.h"
#include "fvscheme.h"

namespace prism {
class CellVisitor {
  public:
    CellVisitor(const mesh::PMesh& mesh, const LinearSystem& system)
        : _mesh(mesh), _system(system) {}

    void visit_cells();
    template <typename T>
    void add_scheme(T& scheme) {
        _schemes.push_back(std::make_shared<T>(scheme));
    }

  private:
    std::vector<std::shared_ptr<FVScheme>> _schemes;
    const mesh::PMesh& _mesh;
    const LinearSystem& _system;
};
} // namespace prism
