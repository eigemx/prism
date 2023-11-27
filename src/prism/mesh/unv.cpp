#include "unv.h"

#include <fmt/format.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string_view>

#include "boundary.h"
#include "prism/exceptions.h"
#include "unvpp/unvpp.h"

// Vector of 6 quad faces
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline hex_cell_faces(const std::vector<std::size_t>& c) -> T {
    T hfaces(6);

    hfaces[0] = {c[1], c[2], c[6], c[5]};
    hfaces[1] = {c[0], c[4], c[7], c[3]};
    hfaces[2] = {c[2], c[3], c[7], c[6]};
    hfaces[3] = {c[0], c[1], c[5], c[4]};
    hfaces[4] = {c[4], c[5], c[6], c[7]};
    hfaces[5] = {c[0], c[3], c[2], c[1]};

    return hfaces;
}

// Vector of 4 triangular faces
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline tetra_cell_faces(const std::vector<std::size_t>& c) -> T {
    T tfaces(4);

    tfaces[0] = {c[0], c[2], c[1]};
    tfaces[1] = {c[1], c[2], c[3]};
    tfaces[2] = {c[0], c[1], c[3]};
    tfaces[3] = {c[0], c[3], c[2]};

    return tfaces;
}

// Vector of 5 faces (3 quads & 2 triangels)
template <typename T = std::vector<std::vector<std::size_t>>>
auto inline wedge_cell_faces(const std::vector<std::size_t>& c) -> T {
    T tfaces(5);

    tfaces[0] = {c[0], c[2], c[1]};
    tfaces[1] = {c[3], c[4], c[5]};
    tfaces[2] = {c[0], c[1], c[4], c[3]};
    tfaces[3] = {c[0], c[3], c[5], c[2]};
    tfaces[4] = {c[1], c[2], c[5], c[4]};

    return tfaces;
}

namespace prism::mesh {
UnvToPMeshConverter::UnvToPMeshConverter(const std::filesystem::path& filename)
    : _filename(filename), unv_mesh(std::make_unique<unvpp::Mesh>(unvpp::read(filename))) {
    // check that mesh has non-zero vertices
    if (unv_mesh->vertices().empty()) {
        throw error::InvalidMesh(fmt::format(
            "UnvToPMeshConverter(): file {} has zero vertices, please check your mesh.",
            filename.c_str()));
    }

    // convert vertices to prism::Vector3d
    for (const auto& v : unv_mesh->vertices()) {
        _vertices.emplace_back(Vector3d(v[0], v[1], v[2]));
    }

    _faces_lookup_trie = std::make_unique<FacesLookupTrie>(unv_mesh->vertices().size());

    process_cells();
    process_groups();

    // clear unv_mesh to save memory
    unv_mesh.release();
}

auto UnvToPMeshConverter::to_pmesh() -> PMesh {
    // TODO: accept custom boundary file path as an argument
    std::vector<std::string_view> boundary_names;

    for (const auto& [name, _] : _boundary_name_to_faces_map) {
        boundary_names.push_back(name);
    }

    // get parent directory from _filename
    auto parent_dir = _filename.parent_path();

    // TODO: accept custom boundary file path as an argument
    auto boundary_patches = read_boundary_file(parent_dir / "boundary.txt", boundary_names);

    // for each boundary patch, set the faces that belong to it
    std::size_t boundary_index {0};
    for (auto& patch : boundary_patches) {
        auto& faces_ids = _boundary_name_to_faces_map.at(patch.name());

        for (auto face_id : faces_ids) {
            _faces[face_id].set_boundary_patch_id(boundary_index);
        }
        ++boundary_index;
    }

    // move all resources to PMesh constructor, UnvToPMeshConverter is no longer needed
    return {std::move(_vertices),
            std::move(_cells),
            std::move(_faces),
            std::move(boundary_patches),
            std::move(_boundary_faces),
            std::move(_interior_faces)};
}

void UnvToPMeshConverter::process_cells() {
    std::size_t unv_element_counter {0};

    for (const auto& element : unv_mesh->elements().value()) {
        switch (element.type()) {
            // we don't care about lines
            case unvpp::ElementType::Line:
                break;

            // boundary faces
            case unvpp::ElementType::Quad:
            case unvpp::ElementType::Triangle: {
                auto boundary_face_id = process_boundary_face(element);

                // we will need UNV id of this boundary face later when we process groups
                _unv_id_to_bface_index_map[unv_element_counter] =
                    // current face has a UNV id equals unv_element_counter, and not yet found in a group
                    std::make_pair(boundary_face_id, false);
                break;
            }

            // cells
            case unvpp::ElementType::Hex:
            case unvpp::ElementType::Tetra:
            case unvpp::ElementType::Wedge:
                process_cell(element);
                break;

            // We shouldn't reach this, because we should've exhausted all possible element types
            // but, you never know!
            default:
                throw error::InvalidMesh(
                    "UnvToPMeshConverter::process_cells(): Input UNV mesh contains an "
                    "unsupported element type!");
                break;
        }

        unv_element_counter++;
    }

    if (_cells.empty()) {
        throw error::InvalidMesh(
            "UnvToPMeshConverter::process_cells(): Input UNV mesh is either empty or two "
            "dimensional. Prism acceptes only 3D meshes");
    }
}

void UnvToPMeshConverter::process_cell(const unvpp::Element& element) {
    // keep track of face ids inside the current cell
    std::vector<std::size_t> cell_faces_ids;

    switch (element.type()) {
        case unvpp::ElementType::Hex: {
            // reserve 6 faces
            cell_faces_ids.reserve(6);
            for (auto& face_vertices_ids : hex_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices_ids);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        case unvpp::ElementType::Tetra: {
            // reserve 4 faces
            cell_faces_ids.reserve(4);
            for (auto& face_vertices_ids : tetra_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices_ids);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        case unvpp::ElementType::Wedge: {
            // reserve 5 faces
            cell_faces_ids.reserve(5);
            for (auto& face_vertices_ids : wedge_cell_faces(element.vertices_ids())) {
                auto face_id = process_face(face_vertices_ids);
                cell_faces_ids.push_back(face_id);
            }
            break;
        }

        default:
            break;
    }

    _cells.emplace_back(
        _faces, std::move(cell_faces_ids), element.vertices_ids(), _cell_id_counter);
    _cell_id_counter++;
}

auto UnvToPMeshConverter::process_face(std::vector<std::size_t>& face_vertices_ids)
    -> std::size_t {
    // sort the face indices to be used in face lookup trie
    auto sorted_face_vertices_ids = face_vertices_ids;
    std::sort(sorted_face_vertices_ids.begin(), sorted_face_vertices_ids.end());

    // lookup the current face, to check if it has been processed before
    auto face_id_opt {face_index(sorted_face_vertices_ids)};

    if (face_id_opt.has_value()) {
        // face exists
        auto face_id {face_id_opt.value()};
        auto& face {_faces[face_id]};

        if (face.has_owner()) {
            // face has been visited before and is owned
            // we update its neighbor and return
            face.set_neighbor(_cell_id_counter);
        } else {
            // face is not owned but has been processed, then it's a boundary face
            // it has been seen already by process_boundary_face()
            // we set the current cell as its owner and return
            face.set_owner(_cell_id_counter);
            _boundary_faces.push_back(face_id);
        }
        return face_id;
    }

    // this a new interior face, we need to call Face() constructor
    // and push the new face to `_faces` vector
    auto face_id {_face_id_counter++};

    std::vector<Vector3d> face_vertices;

    for (const auto& vertex_id : face_vertices_ids) {
        const auto& vec = _vertices[vertex_id];
        face_vertices.push_back(vec);
    }

    auto new_face {Face(face_vertices, std::move(face_vertices_ids))};

    new_face.set_owner(_cell_id_counter);
    new_face.id() = face_id;

    _faces.emplace_back(std::move(new_face));

    _faces_lookup_trie->insert(sorted_face_vertices_ids, face_id);
    _interior_faces.push_back(face_id);

    return face_id;
}

auto UnvToPMeshConverter::process_boundary_face(const unvpp::Element& boundary_face)
    -> std::size_t {
    // sort the face indices to be used as a std::map key for searching
    auto sorted_face_vertices = boundary_face.vertices_ids();
    std::sort(sorted_face_vertices.begin(), sorted_face_vertices.end());

    // this a new face, we need to call Face() constructor
    // and push the new face to `_faces` vector
    auto face_id = _face_id_counter++;

    std::vector<Vector3d> face_vertices;

    for (const auto& vertex_id : boundary_face.vertices_ids()) {
        const auto& vec = _vertices[vertex_id];
        face_vertices.push_back(vec);
    }

    auto new_face = Face(face_vertices, boundary_face.vertices_ids());

    // this face does not yet have an owner, this will be set later by process_face()
    new_face.id() = face_id;
    _faces.emplace_back(std::move(new_face));
    _faces_lookup_trie->insert(sorted_face_vertices, face_id);

    return face_id;
}

auto UnvToPMeshConverter::face_index(const std::vector<std::size_t>& face_vertices) const
    -> std::optional<std::size_t> {
    auto res = _faces_lookup_trie->find(face_vertices);
    if (res.has_value()) {
        return res.value();
    }
    return std::nullopt;
}

void UnvToPMeshConverter::process_groups() {
    if (!unv_mesh->groups()) {
        throw error::InvalidMesh(
            "UnvToPMeshConverter::process_groups(): Input UNV mesh has no groups (boundary "
            "patches)!");
    }

    for (const auto& group : unv_mesh->groups().value()) {
        if (group.type() == unvpp::GroupType::Vertex) {
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
            if ((unique_types.find(unvpp::ElementType::Quad) == unique_types.end()) ||
                (unique_types.find(unvpp::ElementType::Triangle) == unique_types.end())) {
                group_dimension_error = true;
            }
        }

        if (group_dimension_error) {
            throw error::InvalidMesh(
                "UnvToPMeshConverter::process_groups(): Input UNV contains boundary patches with "
                "non-face type elements. Prism supports only boundary patches with "
                "two-dimensional elements.");
        }

        // process the group
        process_group(group);
    }

    check_boundary_faces();
}

void UnvToPMeshConverter::process_group(const unvpp::Group& group) {
    auto group_faces = std::vector<std::size_t>();
    group_faces.reserve(group.elements_ids().size());

    // copy the group elements ids to group_faces after transforming them using _unv_id_to_bface_index_map
    for (const auto& element_id : group.elements_ids()) {
        auto& face_data = _unv_id_to_bface_index_map[element_id];
        group_faces.push_back(face_data.first);
        face_data.second = true;
    }

    _boundary_name_to_faces_map.insert({group.name(), std::move(group_faces)});
}

void UnvToPMeshConverter::check_boundary_faces() {
    std::size_t undefined_boundary_faces_count = 0;

    for (const auto& [unv_id, face_data] : _unv_id_to_bface_index_map) {
        if (!face_data.second) {
            undefined_boundary_faces_count++;
        }
    }

    if (undefined_boundary_faces_count > 0) {
        throw error::InvalidMesh(fmt::format(
            "UnvToPMeshConverter::check_boundary_faces(): Input UNV mesh has {} "
            "boundary faces that are not part of any boundary patches. Please check your mesh.",
            undefined_boundary_faces_count));
    }
}
} // namespace prism::mesh