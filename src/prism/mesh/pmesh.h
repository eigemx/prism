#pragma once

#include <cmath>

#include "boundary.h"
#include "cell.h"
#include "face.h"

namespace prism::mesh {

class PMesh {
public:
    PMesh() = delete;
    PMesh(std::vector<Cell>&& cells, std::vector<Face>&& faces);
    PMesh(std::vector<Cell>&& cells, std::vector<Face>&& faces,
          BoundaryConditions&& boundary_patches);

    [[nodiscard]] auto inline cells() const -> const std::vector<Cell>& { return _cells; }
    [[nodiscard]] auto inline faces() const -> const std::vector<Face>& { return _faces; }
    [[nodiscard]] auto inline boundary_conditions() const -> const BoundaryConditions& {
        return _boundary_conditions;
    }
    void set_boundary_conditions(BoundaryConditions&& boundary_patches);
    void set_boundary_conditions(const BoundaryConditions& boundary_patches);

    auto non_ortho(std::size_t face_id) const -> double;
    auto non_ortho(const Face& face) const -> double;


private:
    std::vector<Cell> _cells;
    std::vector<Face> _faces;
    BoundaryConditions _boundary_conditions;
    const double pi {std::atan(1) * 4};
};
} // namespace prism::mesh