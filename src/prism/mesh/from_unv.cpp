#include "from_unv.h"

#include <algorithm>

namespace prism::mesh {
UnvToPMesh::UnvToPMesh(const std::string& filename) {
    unv_mesh = unv::read(filename);
    process_cells();
}

void UnvToPMesh::process_cells() {
    for (auto& element : unv_mesh.elements.value()) {
        switch (element.type) {
            case unv::ElementType::Hex:
                process_hex_cell(element);

                break;

            default:
                break;
        }
    }
}

void UnvToPMesh::process_hex_cell(const unv::Element& cell_vertices) {
    // keep track of face ids inside the current hexagonal cell
    std::vector<std::size_t> cell_face_ids;

    // reserve a total of 6 faces
    cell_face_ids.reserve(6);

    for (auto& quad_face_vertices : hex_cell_faces(cell_vertices.vertices_ids)) {
        auto face_id = process_face(quad_face_vertices);
        cell_face_ids.push_back(face_id);
    }

    cells.emplace_back(Cell {faces, cell_face_ids, current_cell_id});
    current_cell_id++;
}

auto UnvToPMesh::process_face(const std::vector<std::size_t>& face_vertices)
    -> std::size_t {
    auto sorted_face_vertices = face_vertices;
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());
    auto face_index_iter = face_index(sorted_face_vertices);

    if (face_index_iter.has_value()) {
        // face has been visited before, all we need to do is update it's neighbor
        auto face_id = face_index_iter.value()->second;
        faces[face_id].set_neighbor(current_cell_id);
        return face_id;
    }
    // this a new face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = current_face_id++;
    faces.emplace_back(Face(face_vertices, unv_mesh.vertices));

    return face_id;
}

auto UnvToPMesh::face_index(const std::vector<std::size_t>& face_vertices)
    -> std::optional<UnvToPMesh::FaceIndexMap::iterator> {
    auto face_itr = faces_map.find(face_vertices);

    if (face_itr == faces_map.end()) {
        return std::nullopt;
    }

    return face_itr;
}

} // namespace prism::mesh