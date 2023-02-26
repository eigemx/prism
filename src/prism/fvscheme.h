#pragma once

#include "mesh/pmesh.h"

namespace prism {
class FVScheme {
  public:
    virtual void apply(const mesh::Cell& cell, const mesh::Face& face) const;
    virtual auto run_once() const -> bool = 0;

  private:
    virtual void apply_interior(const mesh::Cell& cell, const mesh::Face& face) const = 0;
    virtual void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const = 0;
};

} // namespace prism