#pragma once

#include "prism/linear.h"
#include "prism/mesh/face.h"

namespace prism::scheme {

// Base type for all finite volume schemes
class FVScheme {
  public:
    FVScheme() = default;
    FVScheme(const FVScheme&) = default;
    FVScheme(FVScheme&&) = default;
    auto operator=(const FVScheme&) -> FVScheme& = default;
    auto operator=(FVScheme&&) -> FVScheme& = default;
    virtual ~FVScheme() = default;

    // apply the discretization scheme
    virtual void apply() = 0;

    // returns true if the scheme requires correction. The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    virtual auto needsCorrection() const noexcept -> bool = 0;
};

// Base type for FVSchemes that requires contribution to only the right hand side of the
// discretized linear system.
class PartialScheme : public FVScheme, public LinearSystem {
  public:
    PartialScheme(std::size_t n_cells) : LinearSystem(n_cells, false) {}

  private:
    VectorXd _rhs;
};

// Base type for FVSchemes that requires contribution to both sides of the discretized linear
// system
template <typename Field>
class FullScheme : public FVScheme, public LinearSystem {
  public:
    FullScheme(std::size_t n_cells, bool need_matrix = true)
        : LinearSystem(n_cells, need_matrix) {}

    // returns the conserved transport field
    virtual auto field() -> Field = 0;

  private:
    virtual void apply_interior(const mesh::Face& face) = 0;
    virtual void apply_boundary(const mesh::Face& face) = 0;
};

} // namespace prism::scheme