#pragma once

#include <vector>

#include "face.h"
#include "prism/types.h"

namespace prism::mesh {
class Cell {
  public:
    Cell(const std::vector<Face>& faces,
         std::vector<std::size_t> faces_ids,
         std::vector<std::size_t> vertices_ids,
         std::size_t cell_id);

    auto volume() const noexcept -> double { return _volume; }
    auto center() const noexcept -> const Vector3d& { return _center; }
    auto id() const noexcept -> std::size_t { return _id; }
    auto id() noexcept -> std::size_t& { return _id; }

    auto facesIds() const noexcept -> const std::vector<std::size_t>& { return _faces_ids; }

    auto facesIds() noexcept -> std::vector<std::size_t>& { return _faces_ids; }

    auto verticesIds() const noexcept -> const std::vector<std::size_t>& { return _vertices_ids; }

  private:
    std::size_t _id {0};
    double _volume {0.0};
    Vector3d _center {0.0, 0.0, 0.0};

    std::vector<std::size_t> _faces_ids;
    std::vector<std::size_t> _vertices_ids;
};

} // namespace prism::mesh
