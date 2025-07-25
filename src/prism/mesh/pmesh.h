#pragma once

#include <cmath>
#include <vector>

#include "boundary.h"
#include "cell.h"
#include "face.h"
#include "iterators.h"
#include "prism/types.h"

namespace prism::mesh {

/// TODO: when loading a unv mesh with two patches having the same name, we get an exception. Fix.
class PMesh {
  public:
    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces,
          std::vector<BoundaryPatch> boundary_patches,
          std::vector<FieldInfo> field_infos,
          std::vector<std::size_t> boundary_faces_ids,
          std::vector<std::size_t> interior_faces_ids) noexcept;

    auto vertices() const noexcept -> const std::vector<Vector3d>&;

    auto cells() const noexcept -> const std::vector<Cell>&;
    auto cells() noexcept -> std::vector<Cell>&;
    auto cell(std::size_t cell_id) const -> const Cell&;
    auto cell(std::size_t cell_id) noexcept -> Cell&;

    auto faces() const noexcept -> const std::vector<Face>&;
    auto faces() noexcept -> std::vector<Face>&;
    auto face(std::size_t face_id) const -> const Face&;
    auto face(std::size_t face_id) noexcept -> Face&;

    auto boundaryPatches() const noexcept -> const std::vector<BoundaryPatch>&;
    auto boundaryPatch(const Face& face) const noexcept -> const BoundaryPatch&;
    auto boundaryPatch(const std::string& name) const -> const BoundaryPatch&;

    auto faceBoundaryPatch(std::size_t face_id) const -> const BoundaryPatch&;
    auto faceBoundaryPatch(const Face& face) const -> const BoundaryPatch&;

    auto cellCount() const noexcept -> std::size_t;
    auto faceCount() const noexcept -> std::size_t;
    auto boundaryFaceCount() const noexcept -> std::size_t;
    auto nonEmptyboundaryFaceCount() const noexcept -> std::size_t;

    auto cellsVolumeVector() const noexcept -> const VectorXd&;

    auto otherSharingCell(const Cell& c, const Face& f) const -> const Cell&;

    auto interiorFaces() const -> iterators::InteriorFaces;
    auto boundaryFaces() const -> iterators::BoundaryFaces;
    auto nonEmptyBoundaryFaces() const -> iterators::BoundaryFaces;

    auto fieldsInfo() const noexcept -> const std::vector<FieldInfo>&;

  private:
    std::vector<Vector3d> _vertices;
    std::vector<Cell> _cells;
    std::vector<Face> _faces;

    std::vector<std::size_t> _boundary_faces_ids;
    std::vector<std::size_t> _nonempty_boundary_faces_ids;
    std::vector<std::size_t> _interior_faces_ids;

    std::vector<BoundaryPatch> _boundary_patches;
    std::vector<FieldInfo> _field_infos;

    std::size_t _n_cells {0};
    std::size_t _n_faces {0};
    VectorXd _cells_volume;
};

class PMeshPtr {
  public:
    PMeshPtr(const PMesh* ptr);

  private:
    const PMesh* _ptr;
};

class ToPMeshConverter {
  public:
    /** @brief Converts a mesh to a PMesh object.
     * Any type that implements this interface isn't expected to be used after calling to_pmesh().
     * this is because the implementation of to_pmesh() is expected to move the data from the
     * implementing object to the returned PMesh object.
     */
    ToPMeshConverter() = default;
    ToPMeshConverter(const ToPMeshConverter&) = default;
    ToPMeshConverter(ToPMeshConverter&&) = default;
    auto operator=(const ToPMeshConverter&) -> ToPMeshConverter& = default;
    auto operator=(ToPMeshConverter&&) -> ToPMeshConverter& = default;
    virtual ~ToPMeshConverter() = default;

    virtual auto toPMesh() -> SharedPtr<PMesh> = 0;
};

} // namespace prism::mesh