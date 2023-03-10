#pragma once

#include <unvpp/unvpp.h>

#include <array>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <utility> // std::pair

#include "cell.h"
#include "face.h"
#include "pmesh.h"

namespace prism::mesh {

class UnvToPMesh : public ToPMeshConverter {
  public:
    UnvToPMesh(const std::filesystem::path& filename);

    auto to_pmesh() -> PMesh override;

  private:
    // map a face sorted vertex ids to its index in `this->_faces` vector
    using SortedFaceToIndexMap = std::map<std::vector<std::size_t>, std::size_t>;

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
    void process_cell(const unv::Element& element);

    // faces
    auto process_face(std::vector<std::size_t>& face_vertices) -> std::size_t;
    auto process_boundary_face(const unv::Element& boundary_face) -> std::size_t;

    // groups
    void process_groups();
    void process_group(const unv::Group& group);
    void check_boundary_faces();

    auto face_index(const std::vector<std::size_t>& face_vertices)
        -> std::optional<SortedFaceToIndexMap::iterator>;

    // fields
    std::filesystem::path _filename;
    unv::Mesh unv_mesh;

    SortedFaceToIndexMap _face_to_index_map;
    BoundaryNameToFacesMap _boundary_name_to_faces_map;
    UnvIndexToBFaceIndexMap _unv_id_to_bface_index_map;

    std::vector<Vector3d> _vertices;
    std::vector<Face> _faces;
    std::vector<Cell> _cells;

    std::size_t _cell_id_counter {0};
    std::size_t _face_id_counter {0};
};

} // namespace prism::mesh