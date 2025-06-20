#pragma once

#include <prism/types.h>

#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

/// TODO: When a fields json file contains a wrong patch name (versteeg1DAdvection case when Front
/// in T.json was mistakenly named as Rront), we got a std::out_of_range exception. Fix this.

/// TODO: in fields.json, it seems that we are not doing any thing with field "type", check this.

namespace prism::mesh {

/** @brief Enum class for boundary condition value types
 *
 * Each boundary condition can have a scalar or a vector value or Nil.
 * The Nil value is used for symmetry and empty boundary conditions, where the value is not
 * required. For Scalar & Vector types the value is set using BoundaryConditionValue variant
 * defined below.
 *
 */
enum class BoundaryConditionValueKind {
    Nil, // Empty or symmetry
    Scalar,
    Vector
};

/** @brief Boundary condition value variant
 *
 * This variant is used to store the value of a boundary condition. It can be either a scalar or a
 * vector.
 *
 */
using BoundaryConditionValue = std::variant<double, Vector3d>;

/** @brief Field info class
 *
 * This class contains information about a field, such as its name, type, value kind, value, and
 * gradient scheme. After parsing "fields.json" file, we get a vector of FieldInfo objects.
 *
 */
class FieldInfo {
  public:
    FieldInfo(std::string name, std::string type)
        : _name(std::move(name)), _type(std::move(type)) {}

    FieldInfo(std::string name, std::string type, std::string grad_scheme)
        : _name(std::move(name)), _type(std::move(type)), _grad_scheme(std::move(grad_scheme)) {}

    FieldInfo(std::string name,
              std::string type,
              std::string grad_scheme,
              std::vector<double> units)
        : _name(std::move(name)),
          _type(std::move(type)),
          _grad_scheme(std::move(grad_scheme)),
          _units(std::move(units)) {}

    FieldInfo(std::string name,
              std::string type,
              std::optional<std::string> grad_scheme,
              std::optional<std::vector<double>> units)
        : _name(std::move(name)),
          _type(std::move(type)),
          _grad_scheme(std::move(grad_scheme)),
          _units(std::move(units)) {}

    auto name() const noexcept -> const std::string& { return _name; }
    auto type() const noexcept -> const std::string& { return _type; }
    auto gradScheme() const noexcept -> const std::optional<std::string>& { return _grad_scheme; }
    auto units() const noexcept -> const std::optional<std::vector<double>>& { return _units; }

  private:
    std::string _name;
    std::string _type;
    std::optional<std::string> _grad_scheme;
    std::optional<std::vector<double>> _units;
};

/** @brief Boundary condition class
 *
 * This class defines a boundary condition for a field. It contains the type of the boundary
 * condition, the value of the boundary condition, and the value type (scalar or vector).
 *
 */
class BoundaryCondition {
  public:
    /// TODO: _value should be a std::option<BoundaryConditionValue>, to avoid having to construct
    // Nil-valued boundary conditions with double value 0.0;
    BoundaryCondition(BoundaryConditionValueKind type,
                      BoundaryConditionValue value,
                      std::string bc_type_str);

    auto valueKind() const noexcept -> BoundaryConditionValueKind { return _value_kind; }
    auto value() const noexcept -> const BoundaryConditionValue& { return _value; }
    auto kindString() const noexcept -> const std::string& { return _kind_str; }

  private:
    std::string _kind_str;
    BoundaryConditionValueKind _value_kind;
    BoundaryConditionValue _value;
};

/** @brief Boundary patch class
 *
 * This class defines a boundary patch. It should be set-up during mesh processing, and it can be
 * used later to get the boundary condition for a specific field, using the field name via
 * functions get_bc() or get_scalar_bc() or get_vector_bc().
 *
 */
class BoundaryPatch {
  public:
    BoundaryPatch(std::string name,
                  std::map<std::string, BoundaryCondition> field_name_to_bc_map);
    auto inline name() const noexcept -> const std::string& { return _name; }

    // this method returns the boundary condition for a scalar field given its name
    auto getBoundaryCondition(const std::string& field_name) const -> const BoundaryCondition&;

    auto getScalarBoundaryCondition(const std::string& field_name) const -> double;
    auto getVectorBoundaryCondition(const std::string& field_name) const -> Vector3d;

    auto facesIds() const noexcept -> const std::vector<std::size_t>& { return _faces_ids; }
    auto facesIds() noexcept -> std::vector<std::size_t>& { return _faces_ids; }

    auto inline isEmpty() const noexcept -> bool { return _is_empty; }

  private:
    auto getScalarBCSubfield(const std::string& name) const -> double;

    std::string _name;
    std::map<std::string, BoundaryCondition> _field_name_to_bc_map;
    std::vector<std::size_t> _faces_ids;
    bool _is_empty {false};
};

class MeshBoundary {
  public:
    MeshBoundary(const std::filesystem::path& fields_path);
    auto inline fields() const noexcept -> const std::vector<FieldInfo>& { return _fields; }
    auto inline patches() const noexcept -> const std::vector<BoundaryPatch>& {
        return _boundary_patches;
    }

  private:
    std::vector<FieldInfo> _fields;
    std::vector<BoundaryPatch> _boundary_patches;
};

} // namespace prism::mesh