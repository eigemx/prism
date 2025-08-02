#pragma once

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
    /// TODO: this should be initialized with the field itself, and implement field(), because all
    /// other schemes already implements field() this way and we need to avoid repetition.
    IFullScheme(std::size_t n_cells) : LinearSystem(n_cells) {}
    void apply() override;

    // returns the conserved transport field
    virtual auto field() -> Field = 0;

  private:
    virtual void applyInterior(const mesh::Face& face) = 0;
    virtual void applyBoundary() = 0;
};

template <typename Field>
void IFullScheme<Field>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh->
     * The discretized equation is applied using applyInterior() (per face basis) and
     * applyBoundary() functions.
     *
     */
    applyBoundary();

    const auto& interior_faces = this->field().mesh()->interiorFaces();
    std::for_each(interior_faces.begin(), interior_faces.end(), [this](const mesh::Face& face) {
        applyInterior(face);
    });

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

} // namespace prism::scheme
