#pragma once

#include <cassert>
#include <cmath>

#include "../types.h"
#include "boundary.h"
#include "cell.h"
#include "face.h"

namespace prism::mesh {

class PMesh {
  public:
    PMesh() = delete;
    PMesh(const PMesh& other) = default;
    PMesh(PMesh&& other) noexcept = default;
    auto operator=(const PMesh& other) -> PMesh& = default;
    auto operator=(PMesh&& other) noexcept -> PMesh& = default;
    ~PMesh() noexcept = default;

    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces) noexcept;

    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces,
          std::vector<BoundaryPatch> boundary_patches) noexcept;

    auto inline vertices() const noexcept -> const std::vector<Vector3d>& { return _vertices; }

    auto inline cells() const noexcept -> const std::vector<Cell>& { return _cells; }
    auto inline cells() noexcept -> std::vector<Cell>& { return _cells; }
    auto inline cell(std::size_t cell_id) const noexcept -> const Cell& {
        return _cells[cell_id];
    }
    auto inline cell(std::size_t cell_id) noexcept -> Cell& { return _cells[cell_id]; }

    auto inline faces() const noexcept -> const std::vector<Face>& { return _faces; }
    auto inline faces() noexcept -> std::vector<Face>& { return _faces; }
    auto inline face(std::size_t face_id) const noexcept -> const Face& {
        return _faces[face_id];
    }
    auto inline face(std::size_t face_id) noexcept -> Face& { return _faces[face_id]; }

    auto inline boundary_patches() const noexcept -> const std::vector<BoundaryPatch>& {
        return _boundary_patches;
    }

    auto inline boundary_patch(const Face& face) const noexcept -> const BoundaryPatch& {
        return _boundary_patches[face.boundary_patch_id().value()];
    }

    void set_boundary_conditions(std::vector<BoundaryPatch> boundary_patches);
    void set_boundary_conditions(const std::vector<BoundaryPatch>& boundary_patches);

    auto face_non_ortho(std::size_t face_id) const -> double;
    auto face_non_ortho(const Face& face) const -> double;

    auto face_boundary_patch(std::size_t face_id) const -> const BoundaryPatch&;
    auto face_boundary_patch(const Face& face) const -> const BoundaryPatch&;

    auto n_cells() const noexcept -> std::size_t { return _n_cells; }

    static auto cells_weighting_factor(const Cell& c, const Cell& n, const Face& f) -> double;

    auto inline neighbor_to(const Cell& c, const Face& f) const -> const Cell& {
        assert(f.has_neighbor() && "PMesh::neighbor_to() called on a boundary face!");
        auto n_id = f.owner() == c.id() ? f.neighbor().value() : f.owner();
        return _cells[n_id];
    }

  private:
    std::vector<Vector3d> _vertices;
    std::vector<Cell> _cells;
    std::vector<Face> _faces;
    std::vector<BoundaryPatch> _boundary_patches;
    std::size_t _n_cells {};
};

class ToPMeshConverter {
  public:
    virtual auto to_pmesh() -> PMesh = 0;
};

} // namespace prism::mesh