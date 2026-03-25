#pragma once

#include <prism/types.h>

#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace prism::mesh {

enum class BoundaryConditionValueKind { Nil, Scalar, Vector };

using BoundaryConditionValue = std::variant<double, Vector3d>;

class FieldInfo {
  public:
    FieldInfo(std::string name, std::string type);
    FieldInfo(std::string name, std::string type, std::string grad_scheme);
    FieldInfo(std::string name,
              std::string type,
              std::string grad_scheme,
              std::vector<double> units);
    FieldInfo(std::string name,
              std::string type,
              std::optional<std::string> grad_scheme,
              std::optional<std::vector<double>> units);

    auto name() const noexcept -> const std::string&;
    auto type() const noexcept -> const std::string&;
    auto gradScheme() const noexcept -> const std::optional<std::string>&;
    auto units() const noexcept -> const std::optional<std::vector<double>>&;

  private:
    std::string _name;
    std::string _type;
    std::optional<std::string> _grad_scheme;
    std::optional<std::vector<double>> _units;
};

class BoundaryCondition {
  public:
    BoundaryCondition(BoundaryConditionValueKind type,
                      BoundaryConditionValue value,
                      std::string bc_type_str);

    auto valueKind() const noexcept -> BoundaryConditionValueKind;
    auto value() const noexcept -> const BoundaryConditionValue&;
    auto kindString() const noexcept -> const std::string&;

  private:
    std::string _kind_str;
    BoundaryConditionValueKind _value_kind;
    BoundaryConditionValue _value;
};

class BoundaryPatch {
  public:
    BoundaryPatch(std::string name,
                  std::map<std::string, BoundaryCondition> field_name_to_bc_map);

    auto name() const noexcept -> const std::string&;
    auto getBoundaryCondition(const std::string& field_name) const -> const BoundaryCondition&;
    auto getScalarBoundaryCondition(const std::string& field_name) const -> double;
    auto getVectorBoundaryCondition(const std::string& field_name) const -> Vector3d;
    auto facesIds() const noexcept -> const std::vector<std::size_t>&;
    auto facesIds() noexcept -> std::vector<std::size_t>&;
    auto isEmpty() const noexcept -> bool;

  private:
    auto getScalarBCSubfield(const std::string& name) const -> double;

    std::string _name;
    std::map<std::string, BoundaryCondition> _field_name_to_bc_map;
    std::vector<std::size_t> _faces_ids;
    bool _is_empty {false};
};

class MeshBoundary {
  public:
    explicit MeshBoundary(const std::filesystem::path& fields_path);

    auto fields() const noexcept -> const std::vector<FieldInfo>&;
    auto patches() const noexcept -> const std::vector<BoundaryPatch>&;

  private:
    std::vector<FieldInfo> _fields;
    std::vector<BoundaryPatch> _boundary_patches;
};

} // namespace prism::mesh