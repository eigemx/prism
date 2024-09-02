#pragma once

#include <prism/mesh/pmesh.h>

struct PitzDailyFields {
    std::vector<prism::Vector3d> velocity;
    std::vector<double> pressure;
    std::vector<double> temperature;
};

auto readFields() -> PitzDailyFields;

class FoamMeshToPMeshConverter : public prism::mesh::ToPMeshConverter {
  public:
    auto toPMesh() -> prism::mesh::PMesh override;
};