#pragma once

#include <fmt/color.h>
#include <fmt/core.h>

#include "pmesh.h"

namespace prism::mesh {

auto non_ortho(std::size_t face_a_id, std::size_t face_b_id) -> double;
class PMeshInspector {
public:
    PMeshInspector() = delete;
    PMeshInspector(const PMesh& pmesh) : pmesh(pmesh) {}

private:
    const PMesh& pmesh;
};

} // namespace prism::mesh