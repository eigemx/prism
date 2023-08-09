#pragma once

#include <optional>
#include <utility>

#include "field.h"
#include "linear.h"
#include "mesh/pmesh.h"
#include "print.h"
#include "types.h"

namespace prism {

class FVScheme : public LinearSystem {
  public:
    FVScheme(std::size_t n_cells) : LinearSystem(n_cells) {}

    // apply the discretization scheme
    virtual void apply() = 0;

    // returns true if the scheme requires correction. The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    virtual auto requires_correction() const -> bool { return true; }

    virtual auto field() -> ScalarField& = 0;

  private:
    virtual void apply_interior(const mesh::Face& face) = 0;
    virtual void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) = 0;
};

} // namespace prism