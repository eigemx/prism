#pragma once

#include <filesystem>
#include <string>
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
    FixedGradient,
    Unknown, // for error handling
};

enum class BoundaryConditionValueType {
    Nil, // Empty or symmetry
    Scalar,
    Vector
};

using BoundaryConditionData = std::variant<double, Vector3d>;

class BoundaryCondition {
  public:
    BoundaryCondition(BoundaryConditionValueType type,
                      BoundaryConditionData data,
                      BoundaryPatchType patch_type)
        : _type(type), _data(std::move(data)), _patch_type(patch_type) {}

    auto type() const noexcept -> BoundaryConditionValueType { return _type; }
    auto data() const noexcept -> const BoundaryConditionData& { return _data; }
    auto patch_type() const noexcept -> BoundaryPatchType { return _patch_type; }

  private:
    std::string _name;
    BoundaryConditionValueType _type;
    BoundaryConditionData _data;
    BoundaryPatchType _patch_type;
};

class BoundaryPatch {
  public:
    BoundaryPatch() = delete;

    BoundaryPatch(std::string name, std::map<std::string, BoundaryCondition> field_name_to_bc_map)
        : _name(std::move(name)), _field_name_to_bc_map(std::move(field_name_to_bc_map)) {}

    auto name() const noexcept -> const std::string& { return _name; }

    // this method returns the boundary condition for a scalar field given its name
    auto get_bc(const std::string& field_name) const -> const BoundaryCondition&;

    auto get_scalar_bc(const std::string& field_name) const -> double;
    auto get_vector_bc(const std::string& field_name) const -> Vector3d;

  private:
    std::string _name;
    std::map<std::string, BoundaryCondition> _field_name_to_bc_map;
};

auto read_boundary_data_file(const std::filesystem::path& path,
                             const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch>;

} // namespace prism::mesh