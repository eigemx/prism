#pragma once

#include <string>

#include "prism/boundary.h"


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

} // namespace prism::scheme::boundary