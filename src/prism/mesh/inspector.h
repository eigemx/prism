#pragma once

#include <fmt/color.h>
#include <fmt/core.h>

#include "pmesh.h"

namespace prism::mesh {
class PMeshInspector {
public:
    PMeshInspector() = delete;
    PMeshInspector(const PMesh& pmesh) : pmesh(pmesh) {}

    void report_mesh_stats() const;
    void report_mesh_connectivity() const;
    void report_boundary_patches() const;

private:
    const PMesh& pmesh;
};

} // namespace prism::mesh