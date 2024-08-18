#pragma once

#include <fmt/format.h>

#include <exception>
#include <string>

namespace prism::error {
class NonImplementedBoundaryCondition : public std::exception {
  public:
    NonImplementedBoundaryCondition(const std::string& func_name,
                                    const std::string& patch_name,
                                    const std::string& boundary_condition_type);
    auto what() const noexcept -> const char* override;

  private:
    std::string _message;
};

class InvalidMesh : public std::exception {
  public:
    InvalidMesh(std::string message);
    auto what() const noexcept -> const char* override;

  private:
    std::string _message;
};

} // namespace prism::error
