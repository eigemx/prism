#include "from_unv.h"

#include <algorithm> // for std::sort
#include <stdexcept>

namespace prism::mesh {
UnvToPMesh::UnvToPMesh(const std::string& filename) {
    unv_mesh = unv::read(filename);
    process_cells();

    auto rand_face = faces[100];
    fmt::print("face id: {}, face owner: {}, face neighbor: {}\n", rand_face.id(),
               rand_face.owner(), rand_face.neighbor().value_or(9999));
}

void UnvToPMesh::process_cells() {
    for (auto& element : unv_mesh.elements.value()) {
        // look for 3-dimensional element types
        switch (element.type) {
            // two dimensional elements, nothing to do
            case unv::ElementType::Line:
            case unv::ElementType::Quad:
            case unv::ElementType::Triangle:
                break;

            // Hexagon
            case unv::ElementType::Hex:
                process_hex_cell(element);
                break;

            // Tetrahedron
            case unv::ElementType::Tetra:
                process_tetra_cell(element);
                break;

            // Wedge (prism)
            case unv::ElementType::Wedge:
                process_wedge_cell(element);
                break;

            // We shouldn't reach this, as unvpp lib won't allow it.
            default:
                throw std::runtime_error("Mesh has unsupported UNV element!");
                break;
        }
    }

    if (cells.empty()) {
        throw std::runtime_error(
            "Input UNV mesh is either empty or two dimensional. Prism acceptes only 3D meshes");
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

void UnvToPMesh::process_tetra_cell(const unv::Element& cell_vertices) {
    // keep track of face ids inside the current tetrahedron cell
    std::vector<std::size_t> cell_face_ids;

    // reserve a total of 4 faces
    cell_face_ids.reserve(4);

    for (auto& tri_face_vertices : tetra_cell_faces(cell_vertices.vertices_ids)) {
        auto face_id = process_face(tri_face_vertices);
        cell_face_ids.push_back(face_id);
    }

    cells.emplace_back(Cell {faces, cell_face_ids, current_cell_id});
    current_cell_id++;
}

void UnvToPMesh::process_wedge_cell(const unv::Element& cell_vertices) {
    // keep track of face ids inside the current wedge cell
    std::vector<std::size_t> cell_face_ids;

    // reserve a total of 4 faces
    cell_face_ids.reserve(5);

    for (auto& tri_face_vertices : tetra_cell_faces(cell_vertices.vertices_ids)) {
        auto face_id = process_face(tri_face_vertices);
        cell_face_ids.push_back(face_id);
    }

    cells.emplace_back(Cell {faces, cell_face_ids, current_cell_id});
    current_cell_id++;
}

auto UnvToPMesh::process_face(const std::vector<std::size_t>& face_vertices) -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = face_vertices;
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // lookup the current face, to check if it has been processed before
    auto face_id_iter = face_index(sorted_face_vertices);

    if (face_id_iter.has_value()) {
        // face has been visited before, all we need to do is update its neighbor
        auto face_id = face_id_iter.value()->second;
        faces[face_id].set_neighbor(current_cell_id);
        return face_id;
    }

    // this a new face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = current_face_id++;
    auto new_face = Face(face_vertices, unv_mesh.vertices);
    new_face.set_owner(current_cell_id);
    new_face.set_id(face_id);
    faces.emplace_back(new_face);

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