#pragma once

#include <fmt/format.h>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

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

template <typename Scheme, typename BaseHandler>
class BoundaryHandlersManager {
  public:
    auto get_handler(const std::string& bc) const -> std::shared_ptr<BaseHandler>;

    template <typename Handler>
    void add_handler();

    void remove_handler(const std::string& bc);

  private:
    using InstanceCreatorPtr = std::shared_ptr<IBoundaryHandler> (*)();
    void add_handler(const std::string& bc, InstanceCreatorPtr creator);
    std::map<std::string, InstanceCreatorPtr> _bc_map;
};

template <typename T>
auto create_handler_instance() -> std::shared_ptr<IBoundaryHandler> {
    return std::make_shared<T>();
}

template <typename Scheme, typename BaseHandler>
auto BoundaryHandlersManager<Scheme, BaseHandler>::get_handler(const std::string& bc) const
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

template <typename Scheme, typename BaseHandler>
void BoundaryHandlersManager<Scheme, BaseHandler>::add_handler(const std::string& bc,
                                                               InstanceCreatorPtr creator) {
    if (creator == nullptr) {
        throw std::runtime_error(
            "BoundaryHandlersManager::add_handler() was given a pointer to boundary handler "
            "instance creator function that is null.");
    }
    _bc_map.insert({bc, creator});
}

template <typename Scheme, typename BaseHandler>
template <typename DerivedHandler>
void BoundaryHandlersManager<Scheme, BaseHandler>::add_handler() {
    DerivedHandler temp;
    add_handler(temp.name(), &boundary::create_handler_instance<DerivedHandler>);
}

template <typename Scheme, typename BaseHandler>
void BoundaryHandlersManager<Scheme, BaseHandler>::remove_handler(const std::string& bc) {
    const auto it = _bc_map.find(bc);
    if (it != _bc_map.end()) {
        _bc_map.erase(it);
        return;
    }
    throw std::runtime_error(
        fmt::format("BoundaryHandlersManager::remove_hander() was given a non-existent boundary "
                    "handler name ({})",
                    bc));
}

} // namespace prism::boundary