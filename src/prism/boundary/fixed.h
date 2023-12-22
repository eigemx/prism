#pragma once

#include <string>

#include "boundary.h"
#include "prism/mesh/boundary.h"

namespace prism::boundary {

template <typename Scheme, typename Field>
class Fixed : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Scheme, typename Field>
void Fixed<Scheme, Field>::apply(Scheme& scheme,
                                 const Field& field,
                                 const mesh::BoundaryPatch& patch) const {}

} // namespace prism::boundary