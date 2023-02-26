#pragma once

#include <cmath>

#include "boundary.h"
#include "cell.h"
#include "face.h"

namespace prism::mesh {

class PMesh {
  public:
    PMesh() = delete;
    PMesh(std::vector<Cell> cells, std::vector<Face> faces);
    PMesh(std::vector<Cell> cells, std::vector<Face> faces, BoundaryPatches boundary_patches);

    [[nodiscard]] auto inline cells() const -> const std::vector<Cell>& { return _cells; }
    [[nodiscard]] auto inline faces() const -> const std::vector<Face>& { return _faces; }
    [[nodiscard]] auto inline boundary_conditions() const -> const BoundaryPatches& {
        return _boundary_patches;
    }
    void set_boundary_conditions(BoundaryPatches boundary_patches);
    void set_boundary_conditions(const BoundaryPatches& boundary_patches);

    auto face_non_ortho(std::size_t face_id) const -> double;
    auto face_non_ortho(const Face& face) const -> double;

    auto face_boundary_patch(std::size_t face_id) const -> const BoundaryPatch&;
    auto face_boundary_patch(const Face& face) const -> const BoundaryPatch&;

  private:
    std::vector<Cell> _cells;
    std::vector<Face> _faces;
    BoundaryPatches _boundary_patches;
    const double pi {std::atan(1) * 4};
};

} // namespace prism::mesh