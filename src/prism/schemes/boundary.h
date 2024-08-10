#pragma once

#include <string>

#include "prism/boundary.h"
#include "prism/exceptions.h"
#include "prism/mesh/pmesh.h"
#include "spdlog/spdlog.h"

namespace prism::scheme::boundary {

template <typename Scheme>
class FVSchemeBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) = 0;
};

template <typename Scheme>
class Empty : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override {};
    auto inline name() const -> std::string override { return "empty"; }
};

template <typename Scheme>
class Fixed : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Scheme>
class NoSlip : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <typename Scheme>
class VelocityInlet : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "velocity-intlet"; }
};

template <typename Scheme>
class Outlet : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "outlet"; }
};

template <typename Scheme>
class Symmetry : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename Scheme>
class FixedGradient : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

namespace detail {
template <typename Scheme>
void applyBoundary(const std::string& scheme_name, Scheme& scheme) {
    assert(_scheme.field().has_value());

    auto _phi = scheme.field().value();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundaryPatches()) {
        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi.name());
        auto handler = scheme.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("{}::apply_boundary()", scheme_name), patch.name(), bc.kindString());
        }

        spdlog::debug(
            "{}::apply_boundary(): applying boundary condition type '{}' on "
            "patch '{}'.",
            scheme_name,
            handler->name(),
            patch.name());

        handler->apply(scheme, patch);
    }
}
} // namespace detail

} // namespace prism::scheme::boundary