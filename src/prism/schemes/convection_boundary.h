#pragma once

#include "boundary.h"


namespace prism::scheme::convection {
// forward declarations
template <typename Field>
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <typename F>
class Fixed<convection::IConvection<F>>
    : public FVSchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename F>
class NoSlip<convection::IConvection<F>>
    : public FVSchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <typename F>
class Symmetry<convection::IConvection<F>>
    : public FVSchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename F>
class Outlet<convection::IConvection<F>>
    : public FVSchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

} // namespace prism::scheme::boundary
