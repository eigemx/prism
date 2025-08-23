#pragma once

#include "prism/scheme/boundary.h"

namespace prism::scheme::convection {
// forward declarations
class IAppliedConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <>
class Fixed<convection::IAppliedConvection>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection> {
  public:
    void apply(convection::IAppliedConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <>
class NoSlip<convection::IAppliedConvection>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection> {
  public:
    void apply(convection::IAppliedConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <>
class Symmetry<convection::IAppliedConvection>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection> {
  public:
    void apply(convection::IAppliedConvection& scheme,
               const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <>
class ZeroGradient<convection::IAppliedConvection>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection> {
  public:
    // we treat zero-gradient boundary condition as outlet condition.
    void apply(convection::IAppliedConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

template <>
class Outlet<convection::IAppliedConvection>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection> {
  public:
    void apply(convection::IAppliedConvection& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};


} // namespace prism::scheme::boundary
