#pragma once

#include <fmt/format.h>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "prism/exceptions.h"
#include "prism/mesh/pmesh.h"
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

template <typename BaseHandler>
class BoundaryHandlersManager {
  public:
    auto getHandler(const std::string& bc) const -> std::shared_ptr<BaseHandler>;

    template <typename Handler>
    void addHandler();

    template <typename Handler, typename... Args>
    void addHandler(Args&&... args);

    void removeHandler(const std::string& bc);

  private:
    using InstanceCreatorPtr = std::shared_ptr<IBoundaryHandler> (*)();
    void addHandler(const std::string& bc, InstanceCreatorPtr creator);
    std::map<std::string, InstanceCreatorPtr> _bc_map;
};

template <typename BaseHandler>
class BHManagersProvider {
  public:
    using ManagerType = BoundaryHandlersManager<BaseHandler>;
    virtual auto boundaryHandlersManager() const noexcept -> const ManagerType& {
        return _bh_manager;
    }
    virtual auto boundaryHandlersManager() -> ManagerType& { return _bh_manager; }

  private:
    ManagerType _bh_manager;
};

template <typename T>
auto createHandlerInstance() -> std::shared_ptr<IBoundaryHandler> {
    return std::make_shared<T>();
}

template <typename BaseHandler>
auto BoundaryHandlersManager<BaseHandler>::getHandler(const std::string& bc) const
    -> std::shared_ptr<BaseHandler> {
    if (!_bc_map.contains(bc)) {
        return nullptr;
    }

    auto handler_creator = _bc_map.at(bc); // handler instance creator function
    auto raw_handler = handler_creator();  // std::shared_ptr<IBoundaryHandler>

    // we need to cast std::shared_ptr<IBoundaryHandler> to
    // std::shared_ptr<BaseHandler>
    auto handler = std::dynamic_pointer_cast<BaseHandler>(raw_handler);

    return handler;
}

template <typename BaseHandler>
void BoundaryHandlersManager<BaseHandler>::addHandler(const std::string& bc,
                                                      InstanceCreatorPtr creator) {
    if (creator == nullptr) {
        throw std::runtime_error(
            "BoundaryHandlersManager::addHandler() was given a pointer to boundary handler "
            "instance creator function that is null.");
    }
    _bc_map.insert({bc, creator});
}

template <typename BaseHandler>
template <typename DerivedHandler>
void BoundaryHandlersManager<BaseHandler>::addHandler() {
    DerivedHandler temp;
    addHandler(temp.name(), &boundary::createHandlerInstance<DerivedHandler>);
}

template <typename BaseHandler>
template <typename DerivedHandler, typename... Args>
void BoundaryHandlersManager<BaseHandler>::addHandler(Args&&... args) {
    DerivedHandler temp(std::forward<Args>(args)...);
    addHandler(temp.name(), &boundary::createHandlerInstance<DerivedHandler>);
}
template <typename BaseHandler>
void BoundaryHandlersManager<BaseHandler>::removeHandler(const std::string& bc) {
    const auto it = _bc_map.find(bc);
    if (it != _bc_map.end()) {
        _bc_map.erase(it);
        return;
    }
    throw std::runtime_error(
        fmt::format("BoundaryHandlersManager::removeHandler() was given a non-existent boundary "
                    "handler name ({})",
                    bc));
}

namespace detail {

template <typename Applier>
void applyBoundary(const std::string& applier_name, Applier& applier) {
    auto _phi = applier.field();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundaryPatches()) {
        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi.name());
        auto handler = applier.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("prism::boundary::detail::applyBoundary() applied by {}",
                            applier_name),
                patch.name(),
                bc.kindString());
        }

        spdlog::debug(
            "prism::boundary::detail::applyBoundary(): applying boundary condition type {} on "
            "patch {} applied by {}.",
            handler->name(),
            patch.name(),
            applier_name);

        handler->apply(applier, patch);
    }
}

template <typename Applier>
void applyBoundaryIfExists(const std::string& applier_name, Applier& applier) {
    auto _phi = applier.field();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundaryPatches()) {
        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi.name());
        auto handler = applier.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            spdlog::debug(
                "prism::boundary::detail::applyBoundaryIfExists(): no equation boundary handler "
                "defined for patch {}, applied by {}, "
                "ignoring...",
                patch.name(),
                applier_name);
            continue;
        }

        spdlog::debug(
            "prism::boundary::detail::applyBoundaryIfExists(): applying boundary condition type "
            "{} on "
            "patch {}, applied by {}'.",
            handler->name(),
            patch.name(),
            applier_name);

        handler->apply(applier, patch);
    }
}
} // namespace detail

} // namespace prism::boundary