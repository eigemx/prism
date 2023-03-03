#include "export_vtu.h"

#include <set>
#include <vector>
#include <vtu11/vtu11.hpp>


void export_to_vtu(const prism::mesh::PMesh& pmesh, const std::string& vtu_filename) {
    std::vector<double> points;
    points.reserve(pmesh.vertices().size() * 3);

    for (const auto& vertex : pmesh.vertices()) {
        points.push_back(vertex.x());
        points.push_back(vertex.y());
        points.push_back(vertex.z());
    }

    std::vector<vtu11::VtkIndexType> connectivity;

    std::vector<vtu11::VtkIndexType> offsets;
    offsets.reserve(pmesh.cells().size());

    // set all cell types to 42 (polyhedron)
    std::vector<vtu11::VtkCellType> types(pmesh.cells().size(), 42);

    for (const auto& cell : pmesh.cells()) {
        // add the vertices to the connectivity vector
        connectivity.insert(
            connectivity.end(), cell.vertices_ids().begin(), cell.vertices_ids().end());

        // add the offset
        offsets.push_back(connectivity.size());
    }

    vtu11::Vtu11UnstructuredMesh mesh {points, connectivity, offsets, types};

    // Create some data associated to points and cells
    std::vector<double> pointData(pmesh.vertices().size(), 0.0);
    std::vector<double> cellData(pmesh.cells().size(), 300.0);

    // Create tuples with (name, association, number of components) for each data set
    std::vector<vtu11::DataSetInfo> dataSetInfo {
        {"Temperature", vtu11::DataSetType::PointData, 1},
        {"Conductivity", vtu11::DataSetType::CellData, 1},
    };

    vtu11::writeVtu(vtu_filename, mesh, dataSetInfo, {pointData, cellData}, "Ascii");
}
