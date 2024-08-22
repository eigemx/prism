#pragma once

#include <spdlog/spdlog.h>

#include <cmath>
#include <iterator>

#include "boundary.h"
#include "cell.h"
#include "face.h"
#include "prism/types.h"

namespace prism::mesh {

// TODO: remove const std::vector& members to std::span
namespace detail {
struct FaceIterator {
    using iterator_category = std::input_iterator_tag;
    using value_type = Face;
    using difference_type = std::ptrdiff_t;
    using pointer = Face*;
    using reference = Face&;

    FaceIterator(const std::vector<Face>& faces,
                 const std::vector<std::size_t>& ids,
                 std::size_t position);

    auto operator++() -> FaceIterator&;
    auto operator++(int) -> FaceIterator;
    auto operator*() const -> const Face&;
    auto operator->() const -> const Face*;
    auto operator==(const FaceIterator& other) const -> bool;
    auto operator!=(const FaceIterator& other) const -> bool;

  private:
    const std::vector<Face>& _faces;
    const std::vector<std::size_t>& _ids;
    std::size_t _current {};
};


struct BoundaryFaces {
    BoundaryFaces(const std::vector<Face>& faces,
                  const std::vector<std::size_t>& boundary_faces_ids);

    auto begin() const -> detail::FaceIterator;
    auto end() const -> detail::FaceIterator;

  private:
    const std::vector<Face>& _faces;
    const std::vector<std::size_t>& _boundary_faces_ids;
};

struct InteriorFaces {
    InteriorFaces(const std::vector<Face>& faces,
                  const std::vector<std::size_t>& interior_faces_ids);

    auto begin() const -> detail::FaceIterator;
    auto end() const -> detail::FaceIterator;

  private:
    const std::vector<Face>& _faces;
    const std::vector<std::size_t>& _interior_faces_ids;
};

} // namespace detail

class PMesh {
  public:
    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces,
          std::vector<BoundaryPatch> boundary_patches,
          std::vector<FieldInfo> field_infos,
          std::vector<std::size_t> boundary_faces_ids,
          std::vector<std::size_t> interior_faces_ids) noexcept;

    auto inline vertices() const noexcept -> const std::vector<Vector3d>& { return _vertices; }

    auto inline cells() const noexcept -> const std::vector<Cell>& { return _cells; }
    auto inline cells() noexcept -> std::vector<Cell>& { return _cells; }
    auto inline cell(std::size_t cell_id) const -> const Cell& { return _cells[cell_id]; }
    auto inline cell(std::size_t cell_id) noexcept -> Cell& { return _cells[cell_id]; }

    auto inline faces() const noexcept -> const std::vector<Face>& { return _faces; }
    auto inline faces() noexcept -> std::vector<Face>& { return _faces; }
    auto inline face(std::size_t face_id) const -> const Face& { return _faces[face_id]; }
    auto inline face(std::size_t face_id) noexcept -> Face& { return _faces[face_id]; }

    auto inline boundaryPatches() const noexcept -> const std::vector<BoundaryPatch>& {
        return _boundary_patches;
    }
    auto inline boundaryPatch(const Face& face) const noexcept -> const BoundaryPatch& {
        assert(face.isBoundary() && face.boundaryPatchId().has_value());
        return _boundary_patches[face.boundaryPatchId().value()];
    }

    auto faceBoundaryPatch(std::size_t face_id) const -> const BoundaryPatch&;
    auto faceBoundaryPatch(const Face& face) const -> const BoundaryPatch&;

    auto nCells() const noexcept -> std::size_t { return _n_cells; }
    auto nFaces() const noexcept -> std::size_t { return _n_faces; }

    auto cellsVolumeVector() const noexcept -> const VectorXd& { return _cells_volume; }

    auto otherSharingCell(const Cell& c, const Face& f) const -> const Cell&;

    auto boundaryFaces() const -> detail::BoundaryFaces { return {_faces, _boundary_faces_ids}; }
    auto interiorFaces() const -> detail::InteriorFaces { return {_faces, _interior_faces_ids}; }

    auto fieldsInfo() const noexcept -> const std::vector<FieldInfo>& { return _field_infos; }

  private:
    std::vector<Vector3d> _vertices;
    std::vector<Cell> _cells;
    std::vector<Face> _faces;
    std::vector<std::size_t> _boundary_faces_ids;
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

    virtual auto to_pmesh() -> PMesh = 0;
};

} // namespace prism::mesh