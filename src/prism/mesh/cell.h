#pragma once

#include <vector>

#include "../types.h"
#include "face.h"

namespace prism::mesh {
class Cell {
  public:
    Cell(const std::vector<Face>& faces,
         std::vector<std::size_t> faces_ids,
         std::vector<std::size_t> vertices_ids,
         std::size_t cell_id);

    auto inline volume() const noexcept -> double { return _volume; }
    auto inline center() const noexcept -> const Vector3d& { return _center; }
    auto inline id() const noexcept -> std::size_t { return _id; }

    auto inline faces_ids() const noexcept -> const std::vector<std::size_t>& {
        return _faces_ids;
    }

    auto inline vertices_ids() const noexcept -> const std::vector<std::size_t>& {
        return _vertices_ids;
    }

  private:
    std::size_t _id {0};
    double _volume {0.0};
    Vector3d _center {0.0, 0.0, 0.0};

    std::vector<std::size_t> _faces_ids;
    std::vector<std::size_t> _vertices_ids;
};

} // namespace prism::mesh
