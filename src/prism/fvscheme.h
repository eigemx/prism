#pragma once

#include <optional>
#include <utility>

#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

using CoeffPair = std::pair<std::size_t, double>;

struct AlteredCoeffs {
    double central {0.0};
    std::optional<CoeffPair> neighbor {std::nullopt};
    double b {0.0};
};

class FVScheme {
  public:
    virtual auto apply(const mesh::Cell& cell,
                       const mesh::Face& face,
                       const mesh::PMesh& mesh) const -> AlteredCoeffs {
        if (!face.has_neighbor()) {
            return apply_boundary(cell, face, mesh);
        }
        return apply_interior(cell, face, mesh);
    }

  private:
    virtual auto apply_interior(const mesh::Cell& cell,
                                const mesh::Face& face,
                                const mesh::PMesh& mesh) const -> AlteredCoeffs = 0;
    virtual auto apply_boundary(const mesh::Cell& cell,
                                const mesh::Face& face,
                                const mesh::PMesh& mesh) const -> AlteredCoeffs = 0;
};

} // namespace prism