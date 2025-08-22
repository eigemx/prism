#pragma once

#include "prism/field/ifield.h"
#include "prism/linear.h"
#include "prism/mesh/face.h"

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

    /** @brief Applies discretized diffusion equation to the mesh->
     * The discretized equation is applied using applyInterior() (per face basis) and
     * applyBoundary() functions.
     *
     */
    virtual void apply() = 0;

    // returns true if the scheme requires correction. The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    virtual auto needsCorrection() const noexcept -> bool = 0;
};

template <typename T>
concept ISchemeBased = std::derived_from<T, IScheme>;

// Base type for FVSchemes that requires contribution to only the right hand side of the
// discretized linear system.
class IPartialScheme : public IScheme, public RHSProvider {
  public:
    IPartialScheme(std::size_t n_cells);
};

// Base type for finite volume schemes that requires contribution to both sides of the discretized
// linear system
class IFullScheme : public IScheme, public LinearSystem {
  public:
    IFullScheme(const SharedPtr<field::IScalar>& field);

    void apply() override;

    // returns the conserved transport field
    virtual auto field() -> SharedPtr<field::IScalar>&;

  private:
    virtual void applyInterior(const mesh::Face& face) = 0;
    virtual void applyBoundary() = 0;

    SharedPtr<field::IScalar> _field;
};

} // namespace prism::scheme
