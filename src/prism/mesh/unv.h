#pragma once

#include <unvpp/unvpp.h>

#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility> // std::pair

#include "cell.h"
#include "face.h"
#include "pmesh.h"
#include "prism/mesh/boundary.h"
#include "trie.h"

namespace prism::mesh {

class UnvToPMeshConverter : public ToPMeshConverter {
  public:
    UnvToPMeshConverter(const std::filesystem::path& mesh_path,
                        const std::filesystem::path& boundary_path);

    auto to_pmesh() -> PMesh override;

  private:
    // BoundaryFaceData (we need a better name) is a pair of:
    // 1) the boundary face index in `this->faces()`.
    // 2) whether the face is contained in a defined boundary batch or not.
    // we use the bool value to check if the faces is an "orphan" boundary face (a boundary face
    // that has no boundary patch defined in the input mesh), "true" means that the face is
    // located in a well defined boundary patch, and "false" otherwise.
    using BoundaryFaceData = std::pair<std::size_t, bool>;

    // maps a boundary face Unv element index to BoundaryFaceData defined above
    using UnvIndexToBFaceDataMap = std::map<std::size_t, BoundaryFaceData>;

    // maps a patch name to a set of faces ids in `_faces` vector
    using PatchNameToFacesIdsMap = std::map<std::string, std::vector<std::size_t>>;

    // cells
    void processCells();
    void processCell(const unvpp::Element& element);

    // faces
    auto processFace(std::vector<std::size_t>& face_vertices) -> std::size_t;
    auto processBoundaryFace(const unvpp::Element& boundary_face) -> std::size_t;

    // groups
    void processGroup();
    void processGroup(const unvpp::Group& group);
    void checkBoundaryFaces();

    auto faceIndex(const std::vector<std::size_t>& face_vertices) const
        -> std::optional<std::size_t>;

    // fields
    std::unique_ptr<unvpp::Mesh> _unv_mesh {nullptr};
    PatchNameToFacesIdsMap _boundary_name_to_faces_map;
    UnvIndexToBFaceDataMap _unv_id_to_bface_index_map;
    std::unique_ptr<FacesLookupTrie> _faces_lookup_trie;
    std::vector<Vector3d> _vertices;
    std::vector<Face> _faces;
    std::vector<Cell> _cells;
    std::vector<std::size_t> _boundary_faces;
    std::vector<std::size_t> _interior_faces;
    std::vector<BoundaryPatch> _boundary_patches;
    std::size_t _cell_id_counter {0};
    std::size_t _face_id_counter {0};
};

} // namespace prism::mesh