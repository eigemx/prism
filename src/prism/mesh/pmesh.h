#pragma once

#include <span>
#include <vector>

#include "boundary.h"
#include "cell.h"
#include "face.h"
#include "prism/types.h"

namespace prism::mesh {

/// TODO: when loading a unv mesh with two patches having the same name, we get an exception. Fix.
class PMesh {
  public:
    /** @brief Constructs a PMesh from raw mesh data.
     *
     * The constructor reorders faces so that non-empty boundary patch faces occupy indices
     * [0, N_nonepty), empty boundary patch faces occupy [N_nonempty, N_boundary), and interior
     * faces occupy [N_boundary, N_faces). This matches the OpenFOAM polyMesh convention, making
     * boundary face iteration trivial.
     *
     * @param vertices Mesh vertex coordinates.
     * @param cells Mesh cells.
     * @param faces Mesh faces (will be reordered).
     * @param boundary_patches Boundary patch definitions.
     * @param field_infos Field metadata.
     * @param boundary_faces_ids Indices of boundary faces in the input faces vector.
     * @param interior_faces_ids Indices of interior faces in the input faces vector.
     */
    PMesh(std::vector<Vector3d> vertices,
          std::vector<Cell> cells,
          std::vector<Face> faces,
          std::vector<BoundaryPatch> boundary_patches,
          std::vector<FieldInfo> field_infos,
          std::vector<size_t> boundary_faces_ids,
          std::vector<size_t> interior_faces_ids) noexcept;

    /** @brief Returns the mesh vertices. */
    auto vertices() const noexcept -> const std::vector<Vector3d>&;

    /** @brief Returns the mesh cells (read-only). */
    auto cells() const noexcept -> const std::vector<Cell>&;

    /** @brief Returns the mesh cells (mutable). */
    auto cells() noexcept -> std::vector<Cell>&;

    /** @brief Returns a single cell by its ID.
     * @param cell_id The cell index.
     * @return Reference to the cell. */
    auto cell(size_t cell_id) const -> const Cell&;

    /** @brief Returns a single cell by its ID (mutable).
     * @param cell_id The cell index.
     * @return Mutable reference to the cell. */
    auto cell(size_t cell_id) noexcept -> Cell&;

    /** @brief Returns all faces (read-only).
     *
     * Faces are ordered as: non-empty boundary, empty boundary, interior. */
    auto faces() const noexcept -> const std::vector<Face>&;

    /** @brief Returns all faces (mutable). */
    auto faces() noexcept -> std::vector<Face>&;

    /** @brief Returns a single face by its ID.
     * @param face_id The face index.
     * @return Reference to the face. */
    auto face(size_t face_id) const -> const Face&;

    /** @brief Returns a single face by its ID (mutable).
     * @param face_id The face index.
     * @return Mutable reference to the face. */
    auto face(size_t face_id) noexcept -> Face&;

    /** @brief Returns the boundary patches. */
    auto boundaryPatches() const noexcept -> const std::vector<BoundaryPatch>&;

    /** @brief Returns the boundary patch containing the given face.
     * @param face A boundary face.
     * @return Reference to the owning boundary patch. */
    auto boundaryPatch(const Face& face) const noexcept -> const BoundaryPatch&;

    /** @brief Returns a boundary patch by name.
     * @param name The patch name.
     * @return Reference to the matching patch.
     * @throws std::runtime_error If no patch matches the given name. */
    auto boundaryPatch(const std::string& name) const -> const BoundaryPatch&;

    /** @brief Returns the boundary patch for a given face index.
     * @param face_id Index of a boundary face.
     * @return Reference to the owning boundary patch. */
    auto faceBoundaryPatch(size_t face_id) const -> const BoundaryPatch&;

    /** @brief Returns the boundary patch for a given face.
     * @param face A boundary face.
     * @return Reference to the owning boundary patch. */
    auto faceBoundaryPatch(const Face& face) const -> const BoundaryPatch&;

    /** @brief Returns the number of cells in the mesh. */
    auto cellCount() const noexcept -> size_t;

    /** @brief Returns the total number of faces in the mesh. */
    auto faceCount() const noexcept -> size_t;

    /** @brief Returns the number of boundary faces in the mesh. */
    auto boundaryFaceCount() const noexcept -> size_t;

    /** @brief Returns a pre-computed vector of cell volumes.
     *
     * The vector is indexed by cell ID and contains the volume of each cell. */
    auto cellsVolumeVector() const noexcept -> const VectorXd&;

    /** @brief Returns the neighbouring cell sharing a given interior face.
     *
     * Given a cell and one of its faces, this returns the cell on the other side of the face.
     *
     * @param c A cell.
     * @param f An interior face belonging to cell c.
     * @return The cell sharing face f with cell c. */
    auto otherSharingCell(const Cell& c, const Face& f) const -> const Cell&;

    /** @brief Returns a span over interior faces. */
    auto interiorFaces() const noexcept -> std::span<const Face>;

    /** @brief Returns a span over boundary faces. */
    auto boundaryFaces() const noexcept -> std::span<const Face>;

    /** @brief Returns a span over boundary faces belonging to non-empty patches only. */
    auto nonEmptyBoundaryFaces() const noexcept -> std::span<const Face>;

    /** @brief Returns field metadata for the mesh. */
    auto fieldsInfo() const noexcept -> const std::vector<FieldInfo>&;

  private:
    /** @brief Stores the old-to-new face ID permutation and face counts. */
    struct FacePermutation {
        std::vector<size_t> old_to_new;
        size_t n_boundary;
        size_t n_nonempty;
    };

    /** @brief Builds the old-to-new face ID permutation.
     *
     * Non-empty boundary faces are assigned positions [0, n_nonempty), empty boundary faces
     * are assigned [n_nonempty, n_boundary), and interior faces are assigned [n_boundary, N).
     *
     * @param boundary_faces_ids Indices of boundary faces in the original ordering.
     * @param interior_faces_ids Indices of interior faces in the original ordering.
     * @return The permutation and face counts. */
    auto buildFacePermutation(const std::vector<size_t>& boundary_faces_ids,
                              const std::vector<size_t>& interior_faces_ids) const
        -> FacePermutation;

    /** @brief Reorders _faces according to the given permutation and updates Face::_id.
     * @param old_to_new The permutation mapping old indices to new indices. */
    void reorderFaces(const std::vector<size_t>& old_to_new);

    /** @brief Remaps face IDs in cells and boundary patches using the permutation.
     * @param old_to_new The permutation mapping old indices to new indices. */
    void remapFaceIds(const std::vector<size_t>& old_to_new);

    std::vector<Vector3d> _vertices;
    std::vector<Cell> _cells;
    std::vector<Face> _faces;

    std::vector<BoundaryPatch> _boundary_patches;
    std::vector<FieldInfo> _field_infos;

    size_t _n_cells {0};
    size_t _n_faces {0};
    size_t _n_boundary_faces {0};
    size_t _n_nonempty_boundary_faces {0};
    VectorXd _cells_volume;
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
