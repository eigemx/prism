#pragma once

#include <unvpp/unvpp.h>

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <utility> // std::pair

#include "cell.h"
#include "face.h"
#include "pmesh.h"
#include "trie.h"

namespace prism::mesh {

class UnvToPMeshConverter : public ToPMeshConverter {
  public:
    UnvToPMeshConverter(const std::filesystem::path& filename);

    auto to_pmesh() -> PMesh override;

  private:
    // map a boundary face Unv element index to:
    // 1) its index in `this->faces()`
    // 2) whether the face is contained in a defined boundary batch or not
    // we use the bool value to check if all boundary faces are contained in a
    // defined boundary batch or not.
    using BFaceData = std::pair<std::size_t, bool>;
    using UnvIndexToBFaceIndexMap = std::map<std::size_t, BFaceData>;

    // map a patch name to a set of faces ids in `_faces` vector
    using BoundaryNameToFacesMap = std::map<std::string, std::vector<std::size_t>>;

    // cells
    void process_cells();
    void process_cell(const unvpp::Element& element);

    // faces
    auto process_face(std::vector<std::size_t>& face_vertices) -> std::size_t;
    auto process_boundary_face(const unvpp::Element& boundary_face) -> std::size_t;

    // groups
    void process_groups();
    void process_group(const unvpp::Group& group);
    void check_boundary_faces();

    auto face_index(const std::vector<std::size_t>& face_vertices) const
        -> std::optional<std::size_t>;

    // fields
    std::filesystem::path _filename;
    unvpp::Mesh unv_mesh;

    BoundaryNameToFacesMap _boundary_name_to_faces_map;
    UnvIndexToBFaceIndexMap _unv_id_to_bface_index_map;
    std::unique_ptr<FacesLookupTrie> _faces_lookup_trie;

    std::vector<Vector3d> _vertices;
    std::vector<Face> _faces;
    std::vector<Cell> _cells;
    std::vector<std::size_t> _boundary_faces;
    std::vector<std::size_t> _interior_faces;

    std::size_t _cell_id_counter {0};
    std::size_t _face_id_counter {0};
};

} // namespace prism::mesh