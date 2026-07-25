#pragma once

#include <fmt/format.h>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "prism/exceptions.h"
#include "prism/log.h"
#include "prism/mesh/boundary.h"
#include "prism/types.h"

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
    auto getHandler(const std::string& bc) const -> SharedPtr<BaseHandler>;

    template <typename Handler>
    void addHandler();

    template <typename Handler, typename... Args>
    void addHandler(Args&&... args);

    void removeHandler(const std::string& bc);

  private:
    std::map<std::string, SharedPtr<IBoundaryHandler>> _handler_cache;
};

template <typename BaseHandler>
class BHManagerProvider {
  public:
    using ManagerType = BoundaryHandlersManager<BaseHandler>;
    virtual auto boundaryHandlersManager() const noexcept -> const ManagerType& {
        return _bh_manager;
    }
    virtual auto boundaryHandlersManager() -> ManagerType& { return _bh_manager; }

  private:
    ManagerType _bh_manager;
};

template <typename BaseHandler>
auto BoundaryHandlersManager<BaseHandler>::getHandler(const std::string& bc) const
    -> SharedPtr<BaseHandler> {
    auto it = _handler_cache.find(bc);
    if (it != _handler_cache.end()) {
        return std::dynamic_pointer_cast<BaseHandler>(it->second);
    }
    return nullptr;
}

template <typename BaseHandler>
template <typename DerivedHandler>
void BoundaryHandlersManager<BaseHandler>::addHandler() {
    auto handler = std::make_shared<DerivedHandler>();
    _handler_cache.emplace(handler->name(), std::move(handler));
}

template <typename BaseHandler>
template <typename DerivedHandler, typename... Args>
void BoundaryHandlersManager<BaseHandler>::addHandler(Args&&... args) {
    auto handler = std::make_shared<DerivedHandler>(std::forward<Args>(args)...);
    _handler_cache.emplace(handler->name(), std::move(handler));
}

template <typename BaseHandler>
void BoundaryHandlersManager<BaseHandler>::removeHandler(const std::string& bc) {
    const auto it = _handler_cache.find(bc);
    if (it != _handler_cache.end()) {
        _handler_cache.erase(it);
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
    const auto& mesh = _phi->mesh();

    for (const auto& patch : mesh->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue;
        }

        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi->name());
        auto handler = applier.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("prism::boundary::detail::applyBoundary() applied by {}", applier_name),
                patch.name(),
                bc.kindString());
        }

        log::debug(
            "prism::boundary::detail::applyBoundary(): applying boundary condition type `{}` on "
            "patch `{}` applied by `{}`.",
            handler->name(),
            patch.name(),
            applier_name);

        handler->apply(applier, patch);
    }
}

template <typename Applier>
void applyBoundaryIfExists(const std::string& applier_name, Applier& applier) {
    auto _phi = applier.field();
    const auto& mesh = _phi->mesh();

    for (const auto& patch : mesh->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue;
        }

        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi->name());
        auto handler = applier.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            log::debug(
                "prism::boundary::detail::applyBoundaryIfExists(): no equation boundary handler "
                "defined for patch `{}`, applied by `{}`, "
                "ignoring...",
                patch.name(),
                applier_name);
            continue;
        }

        log::debug(
            "prism::boundary::detail::applyBoundaryIfExists(): applying boundary condition type "
            "`{}` on patch `{}`, applied by `{}`.",
            handler->name(),
            patch.name(),
            applier_name);

        handler->apply(applier, patch);
    }
}
} // namespace detail

} // namespace prism::boundary
