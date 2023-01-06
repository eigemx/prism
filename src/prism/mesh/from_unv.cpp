#include "from_unv.h"

#include <algorithm> // std::sort
#include <stdexcept>

namespace prism::mesh {
UnvToPMesh::UnvToPMesh(const std::string& filename) {
    unv_mesh = unv::read(filename);
    process_cells();
}

void UnvToPMesh::report_mesh_stats() const {
    fmt::print("Mesh stats:\n");
    fmt::print("Number of cells: {}\n", cells.size());
    fmt::print("Number of faces: {}\n", faces.size());
    fmt::print("Number of vertices: {}\n", unv_mesh.vertices.size());

    std::size_t num_boundary_faces {0};
    for (const auto& [face, face_id] : face_to_index_map) {
        if (faces[face_id].neighbor() == std::nullopt) {
            num_boundary_faces++;
        }
    }
    fmt::print("Number of boundary faces: {}\n", num_boundary_faces);
    fmt::print("Number of internal faces: {}\n", faces.size() - num_boundary_faces);
}

void UnvToPMesh::report_mesh_connectivity() const {
    // for each face in the mesh, print who owns the face and who is its neighbor (if any)
    for (const auto& [face, face_id] : face_to_index_map) {
        fmt::print("Face {} is owned by cell {} and its neighbor is cell {}\n", face_id,
                   faces[face_id].owner(), faces[face_id].neighbor().value_or(9999));
    }
}

void UnvToPMesh::process_cells() {
    std::size_t unv_element_counter {0};

    for (auto& element : unv_mesh.elements.value()) {
        switch (element.type()) {
            case unv::ElementType::Line:
                break;

            case unv::ElementType::Quad:
            case unv::ElementType::Triangle: {
                auto boundary_face_id = process_boundary_face(element);

                // we will need UNV id of this boundary face later when we process groups
                face_id_to_unv_index_map[boundary_face_id] = unv_element_counter;
                break;
            }

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
                throw std::runtime_error("Input UNV mesh has an unsupported element type!");
                break;
        }

        unv_element_counter++;
    }

    if (cells.empty()) {
        throw std::runtime_error(
            "Input UNV mesh is either empty or two dimensional. Prism acceptes only 3D meshes");
    }
}

void UnvToPMesh::process_hex_cell(const unv::Element& element) {
    // keep track of face ids inside the current hexagonal cell
    std::vector<std::size_t> cell_faces_ids;

    // reserve a total of 6 faces
    cell_faces_ids.reserve(6);

    for (auto& quad_face_vertices : hex_cell_faces(element.vertices_ids())) {
        auto face_id = process_face(quad_face_vertices);
        cell_faces_ids.push_back(face_id);
    }

    cells.emplace_back(faces, std::move(cell_faces_ids), current_cell_id);
    current_cell_id++;
}

void UnvToPMesh::process_tetra_cell(const unv::Element& element) {
    // keep track of face ids inside the current tetrahedron cell
    std::vector<std::size_t> cell_faces_ids;

    // reserve a total of 4 faces
    cell_faces_ids.reserve(4);

    for (auto& tri_face_vertices : tetra_cell_faces(element.vertices_ids())) {
        auto face_id = process_face(tri_face_vertices);
        cell_faces_ids.push_back(face_id);
    }

    cells.emplace_back(faces, std::move(cell_faces_ids), current_cell_id);
    current_cell_id++;
}

void UnvToPMesh::process_wedge_cell(const unv::Element& element) {
    // keep track of face ids inside the current wedge cell
    std::vector<std::size_t> cell_faces_ids;

    // reserve a total of 4 faces
    cell_faces_ids.reserve(5);

    for (auto& tri_face_vertices : tetra_cell_faces(element.vertices_ids())) {
        auto face_id = process_face(tri_face_vertices);
        cell_faces_ids.push_back(face_id);
    }

    cells.emplace_back(faces, std::move(cell_faces_ids), current_cell_id);
    current_cell_id++;
}

auto UnvToPMesh::process_face(const std::vector<std::size_t>& face_vertices) -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = face_vertices;
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // lookup the current face, to check if it has been processed before
    auto face_id_iter = face_index(sorted_face_vertices);

    if (face_id_iter.has_value()) {
        auto face_id = face_id_iter.value()->second;
        auto& face = faces[face_id];

        if (face.has_owner()) {
            // face has been visited before and is owned
            // we update its neighbor and return
            face.set_neighbor(current_cell_id);

        } else {
            // face is not owned but has been processed, then it's a boundary face
            // we set the current cell as its owner and return
            face.set_owner(current_cell_id);
        }

        return face_id;
    }

    // this a new interior face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = current_face_id++;
    auto new_face = Face(face_vertices, unv_mesh.vertices);

    new_face.set_owner(current_cell_id);
    new_face.set_id(face_id);

    faces.emplace_back(new_face);

    face_to_index_map.insert({sorted_face_vertices, face_id});

    return face_id;
}

auto UnvToPMesh::process_boundary_face(unv::Element& boundary_face) -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = boundary_face.vertices_ids();
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // this a new face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = current_face_id++;
    auto new_face = Face(boundary_face.vertices_ids(), unv_mesh.vertices);

    // this face does not yet has an owner, this will be set later by process_face()
    new_face.set_id(face_id);
    faces.emplace_back(new_face);
    face_to_index_map.insert({sorted_face_vertices, face_id});

    return face_id;
}

auto UnvToPMesh::face_index(const std::vector<std::size_t>& face_vertices)
    -> std::optional<UnvToPMesh::FaceToIndexMap::iterator> {
    auto face_itr = face_to_index_map.find(face_vertices);

    if (face_itr == face_to_index_map.end()) {
        return std::nullopt;
    }

    return face_itr;
}

void UnvToPMesh::process_groups() {
    if (!unv_mesh.groups.has_value()) {
        throw std::runtime_error("Input UNV mesh has no groups (boundary patches)!");
    }

    for (auto& group : unv_mesh.groups.value()) {
        if (group.type() == unv::GroupType::Vertex) {
            continue;
        }
    }
}

} // namespace prism::mesh