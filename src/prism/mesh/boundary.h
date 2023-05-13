#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <variant>
#include <vector>

#include "../types.h"

namespace prism::mesh {

enum class BoundaryPatchType {
    Fixed,
    Inlet,
    Outlet,
    Symmetry,
    Empty,
    Gradient,
    Unknown, // for error handling
};

enum class BoundaryConditionType { Scalar, Vector };

using BoundaryConditionData = std::variant<double, Vector3d>;

class BoundaryCondition {
  public:
    BoundaryCondition(std::string name, BoundaryConditionType type, BoundaryConditionData data)
        : _name(std::move(name)), _type(type), _data(std::move(data)) {}

    auto name() const noexcept -> const std::string& { return _name; }
    auto type() const noexcept -> BoundaryConditionType { return _type; }
    auto data() const noexcept -> const BoundaryConditionData& { return _data; }

  private:
    std::string _name;
    BoundaryConditionType _type;
    BoundaryConditionData _data;
};

using BoundaryConditions = std::vector<BoundaryCondition>;

class BoundaryPatch {
  public:
    BoundaryPatch() = delete;

    BoundaryPatch(std::string name, BoundaryConditions attributes, BoundaryPatchType type)
        : _name(std::move(name)), _bcs(std::move(attributes)), _type(type) {};

    auto name() const noexcept -> const std::string& { return _name; }
    auto type() const noexcept -> BoundaryPatchType { return _type; }
    auto boundary_conditions() const noexcept -> const BoundaryConditions& { return _bcs; }
    auto get_scalar_bc(const std::string& name) const -> double;
    auto get_vector_bc(const std::string& name) const -> Vector3d;

  private:
    std::string _name;
    BoundaryPatchType _type;
    BoundaryConditions _bcs;
};

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch>;

} // namespace prism::mesh