#include "foam_json_parser.h"

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
    auto fields_file = path("tests/cases/pitzDailyFoam/fields_internal_fields.json");
    auto doc = fileToJson(fields_file);

    auto velocity = doc["U"];
    auto temperature = doc["T"].get<std::vector<double>>();
    auto pressure = doc["P"].get<std::vector<double>>();

    std::vector<prism::Vector3d> velocity_vec;
    std::vector<double> pressure_vec;
    std::vector<double> temperature_vec;

    for (const auto& v : velocity) {
        velocity_vec.emplace_back(v[0], v[1], v[2]);
    }

    for (const auto& p : pressure) {
        pressure_vec.emplace_back(p);
    }

    for (const auto& t : temperature) {
        temperature_vec.emplace_back(t);
    }

    return {velocity_vec, pressure_vec, temperature_vec};
}

auto readFoamMesh() -> FoamMesh {
    auto mesh_file = path("tests/cases/pitzDailyFoam/mesh.json");
    auto doc = fileToJson(mesh_file);
    auto boundary = doc["boundary"];

    std::vector<prism::Vector3d> points_vec;
    for (const auto& p : doc["points"]) {
        points_vec.emplace_back(p[0], p[1], p[2]);
    }

    std::vector<std::vector<std::size_t>> faces_vec;
    for (const auto& f : doc["faces"]) {
        std::vector<size_t> face_vec;
        for (const auto& i : f) {
            face_vec.emplace_back(i);
        }
        faces_vec.emplace_back(face_vec);
    }

    std::vector<std::size_t> owner_vec;
    for (const auto& o : doc["owner"]) {
        owner_vec.emplace_back(o);
    }

    std::vector<int> neighbour_vec;
    for (const auto& n : doc["neighbour"]) {
        neighbour_vec.emplace_back(n);
    }

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

    FoamMesh foam_mesh;
    foam_mesh.points = std::move(points_vec);
    foam_mesh.faces = std::move(faces_vec);
    foam_mesh.owner = std::move(owner_vec);
    foam_mesh.neighbour = std::move(neighbour_vec);
    foam_mesh.patches = std::move(foam_patches);

    return foam_mesh;
}

auto constructFaces(const FoamMesh& foam_mesh) -> std::vector<prism::mesh::Face> {
    using prism::mesh::Face;

    std::vector<Face> faces;
    faces.reserve(foam_mesh.faces.size());

    std::size_t face_id = 0;

    for (const auto& face_vertices_ids : foam_mesh.faces) {
        std::vector<prism::Vector3d> face_vertices;
        face_vertices.reserve(face_vertices_ids.size());

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
    }
    return faces;
}

auto constructCells(const FoamMesh& foam_mesh, const std::vector<prism::mesh::Face>& faces)
    -> std::vector<prism::mesh::Cell> {
    using prism::mesh::Cell;
    CellIdToFaces cell_id_to_faces;

    for (std::size_t face_id = 0; face_id < foam_mesh.faces.size(); ++face_id) {
        const auto& face_owner_id = foam_mesh.owner[face_id];
        const auto& face_neighbor_id = foam_mesh.neighbour[face_id];

        if (cell_id_to_faces.contains(face_owner_id)) {
            cell_id_to_faces[face_owner_id].emplace_back(face_id);
        } else {
            cell_id_to_faces.insert({face_owner_id, {face_id}});
        }

        if (face_neighbor_id >= 0) {
            if (cell_id_to_faces.contains(face_neighbor_id)) {
                cell_id_to_faces[face_neighbor_id].emplace_back(face_id);
            } else {
                cell_id_to_faces.insert({static_cast<std::size_t>(face_neighbor_id), {face_id}});
            }
        }
    }

    std::vector<Cell> cells;
    cells.reserve(foam_mesh.owner.size());

    for (const auto& [cell_id, faces_ids] : cell_id_to_faces) {
        std::vector<std::size_t> vertices_ids;
        for (const auto& face_id : faces_ids) {
            const auto& face = faces[face_id];
            vertices_ids.insert(
                vertices_ids.end(), face.verticesIds().begin(), face.verticesIds().end());
        }

        cells.emplace_back(faces, vertices_ids, faces_ids, cell_id);
    }
    return cells;
}

auto FoamMeshToPMeshConverter::toPMesh() -> prism::mesh::PMesh {
    auto foam_mesh = readFoamMesh();
    auto faces = constructFaces(foam_mesh);
    auto cells = constructCells(foam_mesh, faces);

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

    FoamMeshToPMeshConverter conv;
    auto pmesh = conv.toPMesh();

    auto P = prism::field::Pressure("P",
                                    pmesh,
                                    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
                                        fields.pressure.data(), fields.pressure.size()));

    prism::export_field_vtu(P, "P.vtu");
}