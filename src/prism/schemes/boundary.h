#pragma once

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "prism/field.h"
#include "prism/mesh/boundary.h"
#include "spdlog/spdlog.h"

namespace prism::boundary {

class AbstractBoundaryHandler {
  public:
    AbstractBoundaryHandler() = default;
    AbstractBoundaryHandler(AbstractBoundaryHandler&) = default;
    AbstractBoundaryHandler(AbstractBoundaryHandler&&) noexcept = default;
    auto operator=(const AbstractBoundaryHandler&) -> AbstractBoundaryHandler& = default;
    auto operator=(AbstractBoundaryHandler&&) noexcept -> AbstractBoundaryHandler& = default;
    virtual ~AbstractBoundaryHandler() = default;
};

template <typename Scheme>
class FVSchemeBoundaryHandler : public AbstractBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const = 0;
};

template <typename Scheme>
class Empty : public FVSchemeBoundaryHandler<Scheme> {
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override {};
    auto inline name() const -> std::string override { return "empty"; }
};

template <typename Scheme>
class Fixed : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Scheme>
class VelocityInlet : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "velocity-intlet"; }
};

template <typename Scheme>
class Outlet : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "outlet"; }
};

template <typename Scheme>
class Symmetry : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename Scheme>
class FixedGradient : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

template <typename Scheme>
class SlipWall : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "slip-wall"; }
};

template <typename Scheme>
class NonSlipWall : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "nonslip-wall"; }
};

template <typename T>
auto create_handler_instance() -> std::shared_ptr<AbstractBoundaryHandler> {
    return std::make_shared<T>();
}


template <typename Scheme>
class BoundaryHandlersManager {
  public:
    using InstanceCreatorPtr = std::shared_ptr<AbstractBoundaryHandler> (*)();

    auto get_handler(const std::string& bc) -> std::optional<InstanceCreatorPtr>;
    void add_handler(const std::string& bc, InstanceCreatorPtr creator);

  private:
    std::map<std::string, InstanceCreatorPtr> _bc_map;
};

template <typename Scheme>
auto BoundaryHandlersManager<Scheme>::get_handler(const std::string& bc)
    -> std::optional<InstanceCreatorPtr> {
    if (!_bc_map.contains(bc)) {
        spdlog::error(
            "BoundaryManager::get_handler() do not have a handler for boundary condition with "
            "name: '{}'",
            bc);
        return std::nullopt;
    }
    return _bc_map[bc];
}

template <typename Scheme>
void BoundaryHandlersManager<Scheme>::add_handler(const std::string& bc,
                                                  InstanceCreatorPtr creator) {
    if (creator == nullptr) {
        throw std::runtime_error(
            "BoundaryHandlersManager::add_handler() was given a pointer to boundary handler "
            "instance creator function that is null.");
    }
    _bc_map.insert({bc, creator});
}


} // namespace prism::boundary