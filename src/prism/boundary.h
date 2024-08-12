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

template <typename Applier, typename BaseHandler>
class BoundaryHandlersManager {
  public:
    auto getHandler(const std::string& bc) const -> std::shared_ptr<BaseHandler>;

    template <typename Handler>
    void addHandler();

    void removeHandler(const std::string& bc);

  private:
    using InstanceCreatorPtr = std::shared_ptr<IBoundaryHandler> (*)();
    void addHandler(const std::string& bc, InstanceCreatorPtr creator);
    std::map<std::string, InstanceCreatorPtr> _bc_map;
};

template <typename T>
auto createHandlerInstance() -> std::shared_ptr<IBoundaryHandler> {
    return std::make_shared<T>();
}

template <typename Applier, typename BaseHandler>
auto BoundaryHandlersManager<Applier, BaseHandler>::getHandler(const std::string& bc) const
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

template <typename Applier, typename BaseHandler>
void BoundaryHandlersManager<Applier, BaseHandler>::addHandler(const std::string& bc,
                                                               InstanceCreatorPtr creator) {
    if (creator == nullptr) {
        throw std::runtime_error(
            "BoundaryHandlersManager::addHandler() was given a pointer to boundary handler "
            "instance creator function that is null.");
    }
    _bc_map.insert({bc, creator});
}

template <typename Applier, typename BaseHandler>
template <typename DerivedHandler>
void BoundaryHandlersManager<Applier, BaseHandler>::addHandler() {
    DerivedHandler temp;
    addHandler(temp.name(), &boundary::createHandlerInstance<DerivedHandler>);
}

template <typename Applier, typename BaseHandler>
void BoundaryHandlersManager<Applier, BaseHandler>::removeHandler(const std::string& bc) {
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
void applyBoundary(const std::string& scheme_name, Applier& scheme) {
    assert(scheme.field().has_value());

    auto _phi = scheme.field();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundaryPatches()) {
        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi.name());
        auto handler = scheme.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("{}::applyBoundary()", scheme_name), patch.name(), bc.kindString());
        }

        spdlog::debug(
            "{}::applyBoundary(): applying boundary condition type '{}' on "
            "patch '{}'.",
            scheme_name,
            handler->name(),
            patch.name());

        handler->apply(scheme, patch);
    }
}

template <typename Applier>
void applyBoundaryEquation(const std::string& scheme_name, Applier& applier) {
    auto _phi = applier.field();
    const mesh::PMesh& mesh = _phi.mesh();

    for (const auto& patch : mesh.boundaryPatches()) {
        const mesh::BoundaryCondition& bc = patch.getBoundaryCondition(_phi.name());
        auto handler = applier.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            spdlog::debug(
                "{}::applyBoundaryIfExists: no equation boundary handler defined for patch {}, "
                "ignoring...",
                scheme_name,
                patch.name());
            continue;
        }

        spdlog::debug(
            "{}::apply_boundary(): applying boundary condition type '{}' on "
            "patch '{}'.",
            scheme_name,
            handler->name(),
            patch.name());

        handler->apply(applier, patch);
    }
}
} // namespace detail

} // namespace prism::boundary