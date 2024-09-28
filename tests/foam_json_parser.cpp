#include "foam_json_parser.h"

#include <fmt/core.h>
#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <vector>

using CellIdToFaces = std::map<std::size_t, std::vector<std::size_t>>;

using json = nlohmann::json;
using std::filesystem::path;


struct FoamPatch {
    std::string name;
    std::size_t nFaces {0};
    std::size_t startFace {0};
};

struct FoamMesh {
    std::vector<prism::Vector3d> points;
    std::vector<std::vector<std::size_t>> faces;
    std::vector<std::size_t> owner;
    std::vector<int> neighbour;
    std::vector<FoamPatch> patches;
};

auto fileToJson(const path& path) -> json {
    auto file = std::ifstream(path);
    return json::parse(file);
}

auto readFields() -> PitzDailyFields {
    /**
     * Reads the internal fields from the json file "inetrnal_fields.json" in the directory
     * "tests/cases/pitzDailyFoam/" and returns them as a PitzDailyFields struct.
     */
    auto fields_file = path("tests/cases/pitzDailyFoam/internal_fields.json");
    auto doc = fileToJson(fields_file);

    auto velocity = doc["U"];
    auto temperature = doc["T"].get<std::vector<double>>();
    auto pressure = doc["P"].get<std::vector<double>>();

    std::vector<prism::Vector3d> velocity_vec;
    for (const auto& v : velocity) {
        velocity_vec.emplace_back(v[0], v[1], v[2]);
    }

    return {velocity_vec, pressure, temperature};
}

auto readBoundaryPatches() -> std::vector<FoamPatch> {
    /**
     * Reads the boundary patches from the json file "boundary.json" in the directory
     * "tests/cases/pitzDailyFoam/" and returns them as a std::vector<FoamPatch>.
     */
    auto boundary_patches_file = path("tests/cases/pitzDailyFoam/boundary.json");
    auto boundary_doc = fileToJson(boundary_patches_file);
    auto patches = boundary_doc["patches"];

    std::vector<FoamPatch> foam_patches;
    for (const auto& patch : patches) {
        FoamPatch foam_patch;
        foam_patch.name = patch["name"].get<std::string>();
        foam_patch.nFaces = patch["nFaces"].get<std::size_t>();
        foam_patch.startFace = patch["startFace"].get<std::size_t>();
        foam_patches.emplace_back(foam_patch);
    }
    return foam_patches;
}


auto readFoamMesh() -> FoamMesh {
    auto mesh_file = path("tests/cases/pitzDailyFoam/mesh.json");
    auto doc = fileToJson(mesh_file);

    // Read points
    std::vector<prism::Vector3d> points_vec;
    for (const auto& p : doc["points"]) {
        points_vec.emplace_back(p[0], p[1], p[2]);
    }

    // Read faces
    std::vector<std::vector<std::size_t>> faces_vec;
    for (const auto& f : doc["faces"]) {
        std::vector<size_t> face_vec;
        for (const auto& i : f) {
            face_vec.emplace_back(i);
        }
        faces_vec.emplace_back(face_vec);
    }

    // Read owner
    std::vector<std::size_t> owner_vec;
    for (const auto& o : doc["owner"]) {
        owner_vec.emplace_back(o);
    }

    // Read neighbour
    std::vector<int> neighbour_vec;
    for (const auto& n : doc["neighbour"]) {
        neighbour_vec.emplace_back(n);
    }

    FoamMesh foam_mesh;
    foam_mesh.points = std::move(points_vec);
    foam_mesh.faces = std::move(faces_vec);
    foam_mesh.owner = std::move(owner_vec);
    foam_mesh.neighbour = std::move(neighbour_vec);
    foam_mesh.patches = readBoundaryPatches();

    return foam_mesh;
}

auto constructFaces(const FoamMesh& foam_mesh) -> std::vector<prism::mesh::Face> {
    using prism::mesh::Face;

    std::vector<Face> faces;
    faces.reserve(foam_mesh.faces.size());

    std::size_t face_id = 0;

    for (const auto& face_vertices_ids : foam_mesh.faces) {
        std::vector<prism::Vector3d> face_vertices;

        for (const auto& id : face_vertices_ids) {
            const auto& point = foam_mesh.points[id];
            face_vertices.emplace_back(point);
        }

        Face face_obj(face_vertices, face_vertices_ids);
        face_obj.id() = face_id;
        face_obj.setOwner(foam_mesh.owner[face_id]);

        if (foam_mesh.neighbour[face_id] >= 0) {
            face_obj.setNeighbor(foam_mesh.neighbour[face_id]);
        }
        faces.emplace_back(face_obj);
        face_id++;
    }
    return faces;
}

auto getCellVtkVertices(const std::vector<std::size_t>& cell_faces_ids,
                        const std::vector<prism::mesh::Face>& faces,
                        const std::vector<prism::Vector3d>& points) -> std::vector<std::size_t> {
    // We need to get the lower face and the upper face of the cell (in x-plane)
    // once we have the lower face, we need to get the following vertices:
    // vertex 1 at (xmin, ymin, zmin)
    // vertex 2 at (xmax, ymin, zmin)
    // vertex 3 at (xmin, ymin, zmax)
    // vertex 4 at (xmax, ymin, zmax)
    // and do the same for the upper face (vertices 5 to 8)
    return {};
}

auto constructCells(const std::vector<prism::mesh::Face>& faces, const FoamMesh& foam_mesh)
    -> std::vector<prism::mesh::Cell> {
    using prism::mesh::Cell;
    using std::size_t;
    CellIdToFaces cell_id_to_faces;

    for (const auto& face : faces) {
        const auto face_owner = face.owner();
        const auto face_id = face.id();

        // Link the face to its owner
        if (cell_id_to_faces.contains(face_owner)) {
            cell_id_to_faces[face_owner].emplace_back(face_id);
        } else {
            cell_id_to_faces.insert({face_owner, {face_id}});
        }

        // Link the face to its neighbor if it has one
        if (face.neighbor().has_value()) {
            const auto face_neighbor = face.neighbor().value();

            if (cell_id_to_faces.contains(face_neighbor)) {
                cell_id_to_faces[face_neighbor].emplace_back(face_id);
            } else {
                cell_id_to_faces.insert({face_neighbor, {face_id}});
            }
        }
    }

    std::vector<Cell> cells;
    for (const auto& [cell_id, faces_ids] : cell_id_to_faces) {
        // sanity check that each cell has exactly 6 faces (pitzDaily is all hexahedral)
        if (faces_ids.size() != 6) {
            throw std::runtime_error("Each cell must have exactly 6 faces");
        }

        auto vertices_ids = getCellVtkVertices(faces_ids, faces, foam_mesh.points);
        // sanity check that each cell has exactly 8 vertices
        if (vertices_ids.size() != 8) {
            throw std::runtime_error("Each cell must have exactly 8 vertices");
        }
        cells.emplace_back(faces, vertices_ids, faces_ids, cell_id);
    }
    return cells;
}

auto FoamMeshToPMeshConverter::toPMesh() -> prism::mesh::PMesh {
    auto foam_mesh = readFoamMesh();
    auto faces = constructFaces(foam_mesh);
    auto cells = constructCells(faces, foam_mesh);

    prism::mesh::MeshBoundary mesh_boundary("tests/cases/pitzDailyFoam/fields.json");
    auto boundary_patches = mesh_boundary.patches();
    const auto& field_infos = mesh_boundary.fields();

    auto foam_boundary = fileToJson("tests/cases/pitzDailyFoam/boundary.json");
    auto patches = foam_boundary["patches"];

    for (const auto& patch : patches) {
        const auto& name = patch["name"].get<std::string>();
        auto pmesh_patch = std::find_if(
            boundary_patches.begin(), boundary_patches.end(), [&name](const auto& patch) {
                return patch.name() == name;
            });

        const std::size_t n_faces = patch["nFaces"].get<std::size_t>();
        const std::size_t start_face = patch["startFace"].get<std::size_t>();

        for (std::size_t face_id = start_face; face_id < start_face + n_faces; ++face_id) {
            auto& faces_ids = pmesh_patch->facesIds();
            faces_ids.push_back(face_id);
        }
    }

    std::vector<std::size_t> boundary_faces_ids;
    for (const auto& patch : boundary_patches) {
        std::cout << "patch name: " << patch.name() << std::endl;
        for (const auto& face_id : patch.facesIds()) {
            boundary_faces_ids.push_back(face_id);
        }
    }

    std::vector<std::size_t> interior_faces_ids;
    for (const auto& face : faces) {
        if (face.isInterior()) {
            interior_faces_ids.push_back(face.id());
        }
    }

    return {std::move(foam_mesh.points),
            std::move(cells),
            std::move(faces),
            std::move(boundary_patches),
            field_infos,
            std::move(boundary_faces_ids),
            std::move(interior_faces_ids)};
}


auto main() -> int {
    auto fields = readFields();
    auto foam_mesh = readFoamMesh();
    auto faces = constructFaces(foam_mesh);
    auto cells = constructCells(faces, foam_mesh);
    fmt::print("Cells size: {}", cells.size());
    /*
        FoamMeshToPMeshConverter conv;
        auto pmesh = conv.toPMesh();

        auto P = prism::field::Pressure("P",
                                        pmesh,
                                        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
                                            fields.pressure.data(), fields.pressure.size()));

        prism::export_field_vtu(P, "P.vtu");
    */
}
