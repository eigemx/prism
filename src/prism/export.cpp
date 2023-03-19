#include "export.h"

#include <set>
#include <vector>
#include <vtu11/vtu11.hpp>

namespace prism {
void export_field(const ScalarField& field, const std::string& file_name) {
    const auto& pmesh = field.mesh();

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

    // set all cell types to 12 (hexahedron)
    std::vector<vtu11::VtkCellType> types(pmesh.cells().size(), 12);

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
    std::vector<double> cellData(pmesh.cells().size());

    // Fill the data cell data
    for (std::size_t i = 0; i < pmesh.cells().size(); ++i) {
        cellData[i] = field[i];
    }

    // Create tuples with (name, association, number of components) for each data set
    std::vector<vtu11::DataSetInfo> dataSetInfo {
        {"Points", vtu11::DataSetType::PointData, 1},
        {"Temperature", vtu11::DataSetType::CellData, 1},
    };

    vtu11::writeVtu(file_name, mesh, dataSetInfo, {pointData, cellData}, "Ascii");
}

} // namespace prism