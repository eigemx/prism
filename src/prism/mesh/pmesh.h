#pragma once

#include <cassert>
#include <cmath>
#include <iterator>

#include "../types.h"
#include "boundary.h"
#include "cell.h"
#include "face.h"

namespace prism::mesh {

class PMesh {
  public:
    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces,
          std::vector<BoundaryPatch> boundary_patches,
          std::vector<std::size_t> boundary_faces_ids,
          std::vector<std::size_t> interior_faces_ids) noexcept;

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

    auto face_boundary_patch(std::size_t face_id) const -> const BoundaryPatch&;
    auto face_boundary_patch(const Face& face) const -> const BoundaryPatch&;

    // TODO: this is only needed in mesh checkers,
    // I don't think it should be part of the PMesh object
    auto face_non_ortho(std::size_t face_id) const -> double;
    auto face_non_ortho(const Face& face) const -> double;

    auto n_cells() const noexcept -> std::size_t { return _n_cells; }

    /**
     * @brief Returns the cell that shares the face with the given cell.
     * 
     * @param c Cell to find the other sharing cell of.
     * @param f Face that is shared by the two cells.
     * @return const Cell& The other cell that shares the face with the given cell.
     */
    auto inline other_sharing_cell(const Cell& c, const Face& f) const -> const Cell& {
        assert(f.has_neighbor() && "PMesh::other_sharing_cell() called on a boundary face!");
        auto n_id = f.owner() == c.id() ? f.neighbor().value() : f.owner();
        return _cells[n_id];
    }

    struct FaceIterator {
        using iterator_category = std::input_iterator_tag;
        using value_type = Face;
        using difference_type = std::ptrdiff_t;
        using pointer = Face*;
        using reference = Face&;

        FaceIterator(const std::vector<Face>& faces,
                     const std::vector<std::size_t>& ids,
                     std::size_t position)
            : _faces(faces), _ids(ids), _current(position) {}

        auto operator++() -> FaceIterator& {
            _current++;
            return *this;
        }

        auto operator++(int) -> FaceIterator {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        auto operator*() const -> const Face& { return _faces[_ids[_current]]; }
        auto operator->() const -> const Face* { return &_faces[_ids[_current]]; }

        auto operator==(const FaceIterator& other) const -> bool {
            return _current == other._current;
        }

        auto operator!=(const FaceIterator& other) const -> bool { return !(*this == other); }

      private:
        const std::vector<Face>& _faces;
        const std::vector<std::size_t>& _ids;
        std::size_t _current {};
    };

    struct BoundaryFaces {
        BoundaryFaces(const std::vector<Face>& faces,
                      const std::vector<std::size_t>& boundary_faces_ids)
            : _faces(faces), _boundary_faces_ids(boundary_faces_ids) {}

        auto begin() const -> FaceIterator { return {_faces, _boundary_faces_ids, 0}; }
        auto end() const -> FaceIterator {
            return {_faces, _boundary_faces_ids, _boundary_faces_ids.size()};
        }

      private:
        const std::vector<Face>& _faces;
        const std::vector<std::size_t>& _boundary_faces_ids;
    };

    struct InteriorFaces {
        InteriorFaces(const std::vector<Face>& faces,
                      const std::vector<std::size_t>& interior_faces_ids)
            : _faces(faces), _interior_faces_ids(interior_faces_ids) {}

        auto begin() const -> FaceIterator { return {_faces, _interior_faces_ids, 0}; }
        auto end() const -> FaceIterator {
            return {_faces, _interior_faces_ids, _interior_faces_ids.size()};
        }

      private:
        const std::vector<Face>& _faces;
        const std::vector<std::size_t>& _interior_faces_ids;
    };

    auto boundary_faces() const -> BoundaryFaces { return {_faces, _boundary_faces_ids}; }
    auto interior_faces() const -> InteriorFaces { return {_faces, _interior_faces_ids}; }

  private:
    std::vector<Vector3d> _vertices;
    std::vector<Cell> _cells;
    std::vector<Face> _faces;
    std::vector<BoundaryPatch> _boundary_patches;
    std::vector<std::size_t> _boundary_faces_ids;
    std::vector<std::size_t> _interior_faces_ids;
    std::size_t _n_cells {};
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