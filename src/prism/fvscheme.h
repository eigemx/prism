#pragma once

#include "mesh/pmesh.h"

namespace prism {
class FVScheme {
  public:
    virtual void apply(const mesh::Cell& cell, const mesh::Face& face) const;
    virtual auto run_once() const -> bool = 0;
    virtual auto is_first_run_complete() const -> bool { return _first_run_completed; };
    virtual auto set_first_run_complete() -> void { _first_run_completed = true; }

  private:
    virtual void apply_interior(const mesh::Cell& cell, const mesh::Face& face) const = 0;
    virtual void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const = 0;
    bool _first_run_completed {false};
};

} // namespace prism