#pragma once

#include <vector>

#include "../types.h"
#include "face.h"

namespace prism::mesh {
class Cell {
public:
    Cell() = delete;
    Cell(const std::vector<Face>& faces, std::vector<std::size_t>& faces_ids,
         std::size_t cell_id);

    auto inline volume() -> double& { return _volume; }
    auto inline center() -> Vector3d& { return _center; }
    auto inline id() -> std::size_t& { return _id; }

private:
    std::size_t _id {0};
    std::size_t vertices_count {0};
    double _volume {0.0};
    Vector3d _center {0.0, 0.0, 0.0};
    std::vector<std::size_t> faces_ids;
};
} // namespace prism::mesh
