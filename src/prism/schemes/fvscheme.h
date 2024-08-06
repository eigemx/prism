#pragma once

#include <optional>

#include "prism/linear.h"
#include "prism/mesh/face.h"

namespace prism::scheme {

template <typename Field>
class FVScheme : public LinearSystem {
  public:
    FVScheme(std::size_t n_cells, bool need_matrix = true) : LinearSystem(n_cells, need_matrix) {}

    // apply the discretization scheme
    virtual void apply() = 0;

    // returns true if the scheme requires correction. The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    virtual auto needsCorrection() const -> bool { return true; }

    // If the scheme contributes to the transport equation main matrix, then this method shall
    // return the transport ScalarField, if not, such as the cases of a constant source term where
    // there is no contribution to the main matrix, then this method should return a null option
    // (the base class FVScheme implements this as the default case)
    virtual auto field() -> std::optional<Field> { return std::nullopt; }

  private:
    // TODO: remove apply boundary
    virtual void apply_interior(const mesh::Face& face) = 0;
    virtual void apply_boundary(const mesh::Face& face) = 0;
};

} // namespace prism::scheme