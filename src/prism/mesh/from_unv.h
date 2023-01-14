#pragma once

#include <fmt/core.h>
#include <unvpp/unvpp.h>

#include <array>
#include <map>
#include <optional>
#include <string>
#include <utility> // std::pair

#include "cell.h"
#include "face.h"
#include "pmesh.h"

namespace prism::mesh {

class UnvToPMesh {
public:
    UnvToPMesh() = delete;
    UnvToPMesh(const std::string& filename);

    auto to_pmesh() -> PMesh;

    void report_mesh_stats() const;
    void report_mesh_connectivity() const;
    void report_boundary_patches() const;

private:
    // map a face sorted vertex ids to its index in `this->faces` vector
    using SortedFaceToIndexMap = std::map<std::vector<std::size_t>, std::size_t>;

    // map a boundary face Unv element index to:
    // 1) its index in `this->faces()`
    // 2) whether the face is contained in a defined boundary batch or not
    // we use the bool value to check if all boundary faces are contained in a
    // defined boundary batch or not.
    using BFaceData = std::pair<std::size_t, bool>;
    using UnvIndexToBFaceIndexMap = std::map<std::size_t, BFaceData>;

    // map a patch name to a set of faces ids in `faces` vector
    using BoundaryPatchToFacesMap = std::map<std::string, std::vector<std::size_t>>;

    // cells
    void process_cells();
    void process_cell(const unv::Element& element);

    // faces
    auto process_face(const std::vector<std::size_t>& face_vertices) -> std::size_t;
    auto process_boundary_face(unv::Element& boundary_face) -> std::size_t;

    // groups
    void process_groups();
    void process_group(const unv::Group& group);
    void check_boundary_faces();

    auto face_index(const std::vector<std::size_t>& face_vertices)
        -> std::optional<SortedFaceToIndexMap::iterator>;

    // fields
    unv::Mesh unv_mesh;

    SortedFaceToIndexMap face_to_index_map;
    BoundaryPatchToFacesMap boundary_name_to_faces_map;
    UnvIndexToBFaceIndexMap unv_id_to_bface_index_map;

    std::vector<Face> faces;
    std::vector<Cell> cells;

    std::size_t cell_id_counter {0};
    std::size_t face_id_counter {0};
};

// Array of 6 quad faces
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline hex_cell_faces(const std::vector<std::size_t>& c) -> T {
    T hfaces(6);

    hfaces[0] = {c[1], c[2], c[6], c[5]};
    hfaces[1] = {c[0], c[4], c[7], c[3]};
    hfaces[2] = {c[3], c[7], c[6], c[2]};
    hfaces[3] = {c[0], c[1], c[5], c[4]};
    hfaces[4] = {c[4], c[5], c[6], c[7]};
    hfaces[5] = {c[0], c[3], c[2], c[1]};

    return hfaces;
}


// Array of 4 triangular faces
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline tetra_cell_faces(const std::vector<std::size_t>& c) -> T {
    T tfaces(4);

    tfaces[0] = {c[0], c[2], c[1]};
    tfaces[1] = {c[1], c[2], c[3]};
    tfaces[2] = {c[0], c[1], c[3]};
    tfaces[3] = {c[0], c[3], c[2]};

    return tfaces;
}

// Array of 5 faces (3 quads & 2 triangels)
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline wedge_cell_faces(const std::vector<std::size_t>& c) -> T {
    T tfaces(5);

    tfaces[0] = {c[0], c[2], c[1]};
    tfaces[1] = {c[3], c[4], c[5]};
    tfaces[2] = {c[3], c[0], c[1], c[4]};
    tfaces[3] = {c[0], c[3], c[5], c[2]};
    tfaces[4] = {c[1], c[2], c[5], c[4]};

    return tfaces;
}

} // namespace prism::mesh