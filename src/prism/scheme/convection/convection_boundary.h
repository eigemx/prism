#pragma once

#include "prism/scheme/boundary.h"

namespace prism::scheme::convection {
// forward declarations
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <>
class Fixed<convection::IConvection> : public ISchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto name() const -> std::string override { return "fixed"; }
};

template <>
class NoSlip<convection::IConvection> : public ISchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto name() const -> std::string override { return "no-slip"; }
};

template <>
class Symmetry<convection::IConvection> : public ISchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override {}
    auto name() const -> std::string override { return "symmetry"; }
};

template <>
class ZeroGradient<convection::IConvection>
    : public ISchemeBoundaryHandler<convection::IConvection> {
  public:
    // we treat zero-gradient boundary condition as outlet condition.
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto name() const -> std::string override { return "zero-gradient"; }
};

template <>
class Outlet<convection::IConvection> : public ISchemeBoundaryHandler<convection::IConvection> {
  public:
    void apply(convection::IConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};


} // namespace prism::scheme::boundary
