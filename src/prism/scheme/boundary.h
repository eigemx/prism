#pragma once

#include <string>

#include "prism/boundary.h"


namespace prism::scheme::boundary {

template <typename Scheme>
class ISchemeBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) = 0;
};

template <typename Scheme>
class Fixed : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Scheme>
class NoSlip : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <typename Scheme>
class VelocityInlet : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "velocity-intlet"; }
};

template <typename Scheme>
class Outlet : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "outlet"; }
};

template <typename Scheme>
class Symmetry : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename Scheme>
class FixedGradient : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

template <typename Scheme>
class ZeroGradient : public ISchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

} // namespace prism::scheme::boundary