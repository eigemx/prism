#pragma once

#include <string>

#include "prism/field.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/face.h"
#include "prism/nonortho/nonortho.h"
#include "prism/schemes/convection.h"
#include "prism/schemes/diffusion.h"

namespace prism::boundary {

template <typename Scheme, typename Field>
class BoundaryCondition {
  public:
    BoundaryCondition() = default;
    BoundaryCondition(BoundaryCondition&) = default;
    BoundaryCondition(BoundaryCondition&&) noexcept = default;
    auto operator=(const BoundaryCondition&) -> BoundaryCondition& = default;
    auto operator=(BoundaryCondition&&) noexcept -> BoundaryCondition& = default;
    virtual ~BoundaryCondition() = default;

    virtual void apply(Scheme& scheme,
                       const Field& field,
                       const mesh::BoundaryPatch& patch) const = 0;

    virtual auto name() const -> std::string = 0;
};


} // namespace prism::boundary