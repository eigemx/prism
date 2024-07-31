#pragma once

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "prism/exceptions.h"
#include "prism/field.h"
#include "prism/mesh/boundary.h"
#include "spdlog/spdlog.h"

namespace prism::boundary {

class IBoundaryHandler {
  public:
    IBoundaryHandler() = default;
    IBoundaryHandler(IBoundaryHandler&) = default;
    IBoundaryHandler(IBoundaryHandler&&) noexcept = default;
    auto operator=(const IBoundaryHandler&) -> IBoundaryHandler& = default;
    auto operator=(IBoundaryHandler&&) noexcept -> IBoundaryHandler& = default;
    virtual ~IBoundaryHandler() = default;
};

template <typename Scheme>
class FVSchemeBoundaryHandler : public IBoundaryHandler {
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

template <typename Scheme>
class SlipWall : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "slip-wall"; }
};

template <typename Scheme>
class NonSlipWall : public FVSchemeBoundaryHandler<Scheme> {
  public:
    void apply(Scheme& scheme, const mesh::BoundaryPatch& patch) override = 0;
    auto inline name() const -> std::string override { return "nonslip-wall"; }
};

template <typename T>
auto create_handler_instance() -> std::shared_ptr<IBoundaryHandler> {
    return std::make_shared<T>();
}


template <typename Scheme>
class BoundaryHandlersManager {
  public:
    auto get_handler(const std::string& bc) -> std::shared_ptr<FVSchemeBoundaryHandler<Scheme>>;

    template <typename Handler>
    void add_handler();

  private:
    using InstanceCreatorPtr = std::shared_ptr<IBoundaryHandler> (*)();
    void add_handler(const std::string& bc, InstanceCreatorPtr creator);
    std::map<std::string, InstanceCreatorPtr> _bc_map;
};


template <typename Scheme>
auto BoundaryHandlersManager<Scheme>::get_handler(const std::string& bc)
    -> std::shared_ptr<FVSchemeBoundaryHandler<Scheme>> {
    if (!_bc_map.contains(bc)) {
        spdlog::error(
            "BoundaryManager::get_handler() do not have a handler for boundary condition with "
            "name: '{}'",
            bc);
        return nullptr;
    }

    auto handler_creator = _bc_map[bc];   // handler instance creator function
    auto raw_handler = handler_creator(); // std::shared_ptr<IBoundaryHandler>

    // we need to cast std::shared_ptr<IBoundaryHandler> to
    // std::shared_ptr<FVSchemeBoundaryHandler<Scheme>>
    using SchemeBoundaryHandler = boundary::FVSchemeBoundaryHandler<Scheme>;
    auto handler = std::dynamic_pointer_cast<SchemeBoundaryHandler>(raw_handler);

    return handler;
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

template <typename Scheme>
template <typename Handler>
void BoundaryHandlersManager<Scheme>::add_handler() {
    Handler temp;
    add_handler(temp.name(), &boundary::create_handler_instance<Handler>);
}

namespace detail {
template <typename Scheme>
void apply_boundary(const std::string& scheme_name, Scheme& scheme) {
    assert(_scheme.field().has_value() &&
           "diffusion::detail::apply_boundary() was called on a scheme that does not have a "
           "valid field");

    auto _phi = scheme.field().value();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundary_patches()) {
        const mesh::BoundaryCondition& bc = patch.get_bc(_phi.name());

        spdlog::debug(
            "{}::apply_boundary(): assigning a boundary handler for boundary "
            "condition type '{}' in patch '{}'.",
            scheme_name,
            bc.kind_string(),
            patch.name());

        auto handler = scheme.bc_manager().get_handler(bc.kind_string());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("{}::apply_boundary()", scheme_name), patch.name(), bc.kind_string());
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

} // namespace prism::boundary