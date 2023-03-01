#include "export_vtu.h"

#include <vector>
#include <vtu11/vtu11.hpp>

void export_to_vtu(const prism::mesh::PMesh& pmesh, const std::string& vtu_filename) {
    std::vector<vtu11::VtkIndexType> offsets;
    offsets.reserve(pmesh.cells().size());

    for (const auto& cell : pmesh.cells()) {
        }
}