#include "from_unv.h"

#include <algorithm> // std::sort
#include <stdexcept>

namespace prism::mesh {
UnvToPMesh::UnvToPMesh(const std::string& filename) {
    unv_mesh = unv::read(filename);
    process_cells();
    process_groups();
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

    report_mesh_connectivity();
    report_boundary_patches();

    // print mesh volume
    double mesh_volume {0.0};
    for (const auto& cell : cells) {
        mesh_volume += cell.volume();
    }
    fmt::print("Mesh volume = {} L^3\n", mesh_volume);

    // print mesh surface area
    double mesh_surface_area {0.0};
    for (const auto& face : faces) {
        if (face.neighbor() == std::nullopt) {
            mesh_surface_area += face.area();
        }
    }
    fmt::print("Mesh surface area: {} L^2\n", mesh_surface_area);
}

void UnvToPMesh::report_mesh_connectivity() const {
    // print number of faces without owners
    std::size_t num_faces_without_owners {0};
    for (const auto& [face, face_id] : face_to_index_map) {
        if (!faces[face_id].has_owner()) {
            num_faces_without_owners++;
        }
    }
    fmt::print("Number of faces without owners: {}\n", num_faces_without_owners);
}

void UnvToPMesh::report_boundary_patches() const {
    // for each boundary patch, print its name and the count of faces it contains
    for (const auto& [patch_name, faces] : boundary_name_to_faces_map) {
        fmt::print("Boundary patch {} contains {} faces\n", patch_name, faces.size());
    }
}

void UnvToPMesh::process_cells() {
    std::size_t unv_element_counter {0};

    for (auto& element : unv_mesh.elements.value()) {
        switch (element.type()) {
            case unv::ElementType::Line:
                break;

            // boundary faces
            case unv::ElementType::Quad:
            case unv::ElementType::Triangle: {
                auto boundary_face_id = process_boundary_face(element);

                // we will need UNV id of this boundary face later when we process groups
                unv_id_to_bface_index_map[unv_element_counter] =
                    // current face has a Unv id of unv_element_counter, and not yet found in a group
                    std::make_pair(boundary_face_id, false);
                break;
            }

            // cells
            case unv::ElementType::Hex:
            case unv::ElementType::Tetra:
            case unv::ElementType::Wedge:
                process_cell(element);
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

void UnvToPMesh::process_cell(const unv::Element& element) {
    // keep track of face ids inside the current cell
    std::vector<std::size_t> cell_faces_ids;

    switch (element.type()) {
        case unv::ElementType::Hex: {
            // reserve a total of 6 faces
            cell_faces_ids.reserve(6);
            for (auto& face_vertices : hex_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        case unv::ElementType::Tetra: {
            // reserve a total of 4 faces
            cell_faces_ids.reserve(4);
            for (auto& face_vertices : tetra_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        case unv::ElementType::Wedge: {
            // reserve a total of 5 faces
            cell_faces_ids.reserve(5);
            for (auto& face_vertices : wedge_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        default:
            break;
    }

    cells.emplace_back(faces, std::move(cell_faces_ids), cell_id_counter);
    cell_id_counter++;
}

auto UnvToPMesh::process_face(const std::vector<std::size_t>& face_vertices) -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = face_vertices; // copy the face vertices to sort in place
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // lookup the current face, to check if it has been processed before
    auto face_id_iter = face_index(sorted_face_vertices);

    if (face_id_iter.has_value()) {
        auto face_id = face_id_iter.value()->second;
        auto& face = faces[face_id];

        if (face.has_owner()) {
            // face has been visited before and is owned
            // we update its neighbor and return
            face.set_neighbor(cell_id_counter);

        } else {
            // face is not owned but has been processed, then it's a boundary face
            // we set the current cell as its owner and return
            face.set_owner(cell_id_counter);
        }

        return face_id;
    }

    // this a new interior face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = face_id_counter++;
    auto new_face = Face(face_vertices, unv_mesh.vertices);

    new_face.set_owner(cell_id_counter);
    new_face.set_id(face_id);

    faces.push_back(std::move(new_face));

    face_to_index_map.insert({sorted_face_vertices, face_id});

    return face_id;
}

auto UnvToPMesh::process_boundary_face(unv::Element& boundary_face) -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = boundary_face.vertices_ids();
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // this a new face, we need to call Face() constructor
    // and push the new face to `faces` vector
    auto face_id = face_id_counter++;
    auto new_face = Face(boundary_face.vertices_ids(), unv_mesh.vertices);

    // this face does not yet has an owner, this will be set later by process_face()
    new_face.set_id(face_id);
    faces.push_back(std::move(new_face));
    face_to_index_map.insert({sorted_face_vertices, face_id});

    return face_id;
}

auto UnvToPMesh::face_index(const std::vector<std::size_t>& face_vertices)
    -> std::optional<UnvToPMesh::SortedFaceToIndexMap::iterator> {
    auto face_itr = face_to_index_map.find(face_vertices);

    if (face_itr == face_to_index_map.end()) {
        return std::nullopt;
    }

    return face_itr;
}

void UnvToPMesh::process_groups() {
    if (!unv_mesh.groups) {
        throw std::runtime_error("Input UNV mesh has no groups (boundary patches)!");
    }

    for (auto& group : unv_mesh.groups.value()) {
        if (group.type() == unv::GroupType::Vertex) {
            continue;
        }
        const auto& unique_types = group.unique_element_types();

        bool group_dimension_error = false;

        // unique elements count should be lower than or equal 2,
        // and should be quad, triangle or both
        if (unique_types.size() > 2) {
            group_dimension_error = true;
        }

        if (unique_types.size() == 2) {
            if ((unique_types.find(unv::ElementType::Quad) == unique_types.end()) ||
                (unique_types.find(unv::ElementType::Triangle) == unique_types.end())) {
                group_dimension_error = true;
            }
        }

        if (group_dimension_error) {
            throw std::runtime_error(
                "Input UNV contains boundary patches with non-face type elements. "
                "Prism supports only boundary patches with two-dimensional elements.");
        }

        // process the group
        process_group(group);
    }

    check_boundary_faces();
}

void UnvToPMesh::process_group(const unv::Group& group) {
    auto group_faces = std::vector<std::size_t>();
    group_faces.reserve(group.elements_ids().size());

    // copy the group elements ids to group_faces after transforming them using unv_id_to_bface_index_map
    for (const auto& element_id : group.elements_ids()) {
        auto& face_data = unv_id_to_bface_index_map[element_id];
        group_faces.push_back(face_data.first);
        face_data.second = true;
    }

    boundary_name_to_faces_map.insert({group.name(), std::move(group_faces)});
}

void UnvToPMesh::check_boundary_faces() {
    std::size_t undefined_boundary_faces_count = 0;

    for (const auto& [unv_id, face_data] : unv_id_to_bface_index_map) {
        if (!face_data.second) {
            undefined_boundary_faces_count++;
        }
    }

    if (undefined_boundary_faces_count > 0) {
        throw std::runtime_error(fmt::format(
            "Input UNV mesh has {} boundary faces that are not part of any boundary patch. "
            "Please check your mesh.",
            undefined_boundary_faces_count));
    }
}
} // namespace prism::mesh