#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <variant>
#include <vector>

#include "prism/types.h"

namespace prism::mesh {

/** @brief Enum class for boundary condition types
 *
 * This enum class defines the types of boundary conditions that can be applied to a boundary
 * patch. Each boundary patch can have at most one boundary condition type per field. For example,
 * a boundary patch called 'outlet' can have a fixed value boundary condition for the pressure
 * field and a fixed gradient boundary condition for the velocity field.
 *
 * This can be defined in boundary.txt file as:
 *
 * [outlet]
 * [outlet.pressure]
 * type = "fixed"
 * value = 0.0
 * [outlet.velocity]
 * type = "gradient"
 * value = 0.0
 *
 */
enum class BoundaryConditionKind {
    Fixed,
    VelocityInlet,
    Outlet,
    Symmetry,
    Empty,
    FixedGradient,
    SlipWall,
    NonSlipWall,
    Unknown, // for error handling
};

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


/** @brief Boundary condition class
 *
 * This class defines a boundary condition for a field. It contains the type of the boundary
 * condition, the value of the boundary condition, and the value type (scalar or vector).
 *
 */
class BoundaryCondition {
  public:
    BoundaryCondition(BoundaryConditionValueKind type,
                      BoundaryConditionValue value,
                      BoundaryConditionKind patch_type,
                      std::string bc_type_str)
        : _value_kind(type),
          _value(std::move(value)),
          _kind(patch_type),
          _kind_str(std::move(bc_type_str)) {}

    auto inline value_kind() const noexcept -> BoundaryConditionValueKind {
        // Should we check if the value is Nil?
        // this always assumes that value is either scalar or vector.
        return _value_kind;
    }
    auto inline value() const noexcept -> const BoundaryConditionValue& { return _value; }
    auto inline kind() const noexcept -> BoundaryConditionKind { return _kind; }
    auto inline kind_string() const noexcept -> const std::string& { return _kind_str; }

  private:
    std::string _kind_str;
    BoundaryConditionValueKind _value_kind;
    BoundaryConditionValue _value;
    BoundaryConditionKind _kind;
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
    BoundaryPatch() = delete;

    BoundaryPatch(std::string name, std::map<std::string, BoundaryCondition> field_name_to_bc_map)
        : _name(std::move(name)), _field_name_to_bc_map(std::move(field_name_to_bc_map)) {}

    auto inline name() const noexcept -> const std::string& { return _name; }

    // this method returns the boundary condition for a scalar field given its name
    auto get_bc(const std::string& field_name) const -> const BoundaryCondition&;

    auto get_scalar_bc(const std::string& field_name) const -> double;
    auto get_vector_bc(const std::string& field_name) const -> Vector3d;

    auto faces_ids() const noexcept -> const std::vector<std::size_t>& { return _faces_ids; }
    auto faces_ids() noexcept -> std::vector<std::size_t>& { return _faces_ids; }

  private:
    auto get_scalar_bc_subfield(const std::string& name) const -> double;

    std::string _name;
    std::map<std::string, BoundaryCondition> _field_name_to_bc_map;
    std::vector<std::size_t> _faces_ids;
};

/** @brief Read boundary file
 *
 * This function reads the boundary file and returns a vector of boundary patches.
 *
 * @param path Path to the boundary file
 * @param boundary_names Vector of expected boundary names to read from the file
 * @return std::vector<BoundaryPatch> Vector of boundary patches
 */
auto read_boundary_file(const std::filesystem::path& path,
                        const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch>;

} // namespace prism::mesh