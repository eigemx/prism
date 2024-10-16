#include "exceptions.h"

namespace prism::error {
NonImplementedBoundaryCondition::NonImplementedBoundaryCondition(
    const std::string& func_name,
    const std::string& patch_name,
    const std::string& boundary_condition_type) {
    _message = fmt::format(
        "In function `{}`: "
        "Non-implemented boundary condition type `{}` for boundary patch: '{}'",
        func_name,
        boundary_condition_type,
        patch_name);
}

auto NonImplementedBoundaryCondition::what() const noexcept -> const char* {
    return _message.c_str();
}

InvalidMesh::InvalidMesh(std::string message) : _message(std::move(message)) {}

auto InvalidMesh::what() const noexcept -> const char* {
    return _message.c_str();
}

} // namespace prism::error