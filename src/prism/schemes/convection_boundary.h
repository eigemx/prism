#pragma once

#include "boundary.h"


namespace prism::scheme::convection {
// forward declarations
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <>
class Fixed<convection::IConvection> : public FVSchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <>
class Symmetry<convection::IConvection>
    : public FVSchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <>
class Outlet<convection::IConvection> : public FVSchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

} // namespace prism::scheme::boundary
