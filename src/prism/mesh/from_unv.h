#pragma once

#include <fmt/core.h>
#include <unvpp/unvpp.h>

#include <array>
#include <map>
#include <optional>
#include <string>
#include <utility>

#include "cell.h"
#include "face.h"

namespace prism::mesh {

class UnvToPMesh {
public:
    UnvToPMesh() = delete;
    UnvToPMesh(const std::string& filename);

    void quick_report() {
        // in this house, we debug by printing.
        fmt::print("Number of faces = {}, number of cells = {}\n", faces.size(), cells.size());
    }

private:
    using FaceIndexMap = std::map<std::vector<std::size_t>, std::size_t>;
    void process_cells();
    void process_hex_cell(const unv::Element& cell_vertices);
    void process_tetra_cell(const unv::Element& cell_vertices);
    void process_wedge_cell(const unv::Element& cell_vertices);
    auto process_face(const std::vector<std::size_t>& face_vertices) -> std::size_t;
    auto face_index(const std::vector<std::size_t>& face_vertices)
        -> std::optional<FaceIndexMap::iterator>;

    unv::Mesh unv_mesh;
    FaceIndexMap faces_map; // maps a sorted face to its index in `faces` vector
    std::vector<Face> faces;
    std::vector<Cell> cells;

    std::size_t current_cell_id {0};
    std::size_t current_face_id {0};
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