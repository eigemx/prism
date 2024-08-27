#pragma once

#include "prism/linear.h"

namespace prism::scheme {

// Base type for all finite volume schemes
class IScheme {
  public:
    IScheme() = default;
    IScheme(const IScheme&) = default;
    IScheme(IScheme&&) = default;
    auto operator=(const IScheme&) -> IScheme& = default;
    auto operator=(IScheme&&) -> IScheme& = default;
    virtual ~IScheme() = default;

    // apply the discretization scheme
    virtual void apply() = 0;

    // returns true if the scheme requires correction. The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    virtual auto needsCorrection() const noexcept -> bool = 0;
};

// Base type for FVSchemes that requires contribution to only the right hand side of the
// discretized linear system.
class IPartialScheme : public IScheme, public RHSProvider {
  public:
    IPartialScheme(std::size_t n_cells) : RHSProvider(n_cells) {}
};

// Base type for FVSchemes that requires contribution to both sides of the discretized linear
// system
template <typename Field>
class IFullScheme : public IScheme, public LinearSystem {
  public:
    IFullScheme(std::size_t n_cells) : LinearSystem(n_cells) {}

    // returns the conserved transport field
    virtual auto field() -> Field = 0;
};

} // namespace prism::scheme