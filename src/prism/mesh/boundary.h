#pragma once

#include <prism/types.h>

#include <filesystem>
#include <string>
#include <variant>
#include <vector>

namespace prism::mesh {

/**
 * @brief Kind of value stored in a boundary condition.
 */
enum class BoundaryConditionValueKind { Nil, Scalar, Vector };

/**
 * @brief The value of a boundary condition, either a scalar or a vector.
 */
using BoundaryConditionValue = std::variant<double, Vector3d>;

/**
 * @brief Represents metadata for a field (name, type, optional gradient scheme and units).
 */
class FieldInfo {
  public:
    /**
     * @brief Constructs a FieldInfo with just name and type.
     * @param name The name of the field.
     * @param type The type of the field.
     */
    FieldInfo(std::string name, std::string type);

    /**
     * @brief Constructs a FieldInfo with all optional parameters.
     * @param name The name of the field.
     * @param type The type of the field.
     * @param grad_scheme Optional gradient scheme.
     * @param units Optional units.
     */
    FieldInfo(std::string name,
              std::string type,
              Optional<std::string> grad_scheme,
              Optional<std::vector<double>> units);

    /**
     * @brief Returns the name of the field.
     * @return const reference to the field name.
     */
    auto name() const noexcept -> const std::string&;

    /**
     * @brief Returns the type of the field.
     * @return const reference to the field type.
     */
    auto type() const noexcept -> const std::string&;

    /**
     * @brief Returns the optional gradient scheme.
     * @return const reference to the optional gradient scheme.
     */
    auto gradScheme() const noexcept -> const Optional<std::string>&;

    /**
     * @brief Returns the optional units.
     * @return const reference to the optional units.
     */
    auto units() const noexcept -> const Optional<std::vector<double>>&;

  private:
    std::string _name;
    std::string _type;
    Optional<std::string> _grad_scheme;
    Optional<std::vector<double>> _units;
};

/**
 * @brief Represents a boundary condition with a value and a kind string.
 */
class BoundaryCondition {
  public:
    /**
     * @brief Constructs a BoundaryCondition.
     * @param type The value kind (Nil, Scalar, or Vector).
     * @param value The boundary condition value.
     * @param bc_type_str The boundary condition type string (e.g., "fixedValue", "zeroGradient").
     */
    BoundaryCondition(BoundaryConditionValueKind type,
                      BoundaryConditionValue value,
                      std::string bc_type_str);

    /**
     * @brief Returns the value kind of this boundary condition.
     * @return The BoundaryConditionValueKind.
     */
    auto valueKind() const noexcept -> BoundaryConditionValueKind;

    /**
     * @brief Returns the value of this boundary condition.
     * @return const reference to the BoundaryConditionValue.
     */
    auto value() const noexcept -> const BoundaryConditionValue&;

    /**
     * @brief Returns the kind string of this boundary condition.
     * @return const reference to the kind string.
     */
    auto kindString() const noexcept -> const std::string&;

  private:
    std::string _kind_str;
    BoundaryConditionValueKind _value_kind;
    BoundaryConditionValue _value;
};

/**
 * @brief Represents a boundary patch with a name and a map of field names to boundary conditions.
 */
class BoundaryPatch {
  public:
    /**
     * @brief Constructs a BoundaryPatch.
     * @param name The name of the boundary patch.
     * @param field_name_to_bc_map Map of field names to their boundary conditions.
     * @throws std::runtime_error If the field name to BC map is empty.
     */
    BoundaryPatch(std::string name,
                  std::map<std::string, BoundaryCondition> field_name_to_bc_map);

    /**
     * @brief Returns the name of the boundary patch.
     * @return const reference to the patch name.
     */
    auto name() const noexcept -> const std::string&;

    /**
     * @brief Gets the boundary condition for a given field name.
     * @param field_name The name of the field.
     * @return A const reference to the BoundaryCondition.
     * @throws std::runtime_error If the field is not found in the patch.
     */
    auto getBoundaryCondition(const std::string& field_name) const -> const BoundaryCondition&;

    /**
     * @brief Gets the scalar boundary condition value for a given field.
     * @param field_name The name of the field.
     * @return The scalar value.
     * @throws std::runtime_error If the field is not found or is not a scalar.
     */
    auto getScalarBoundaryCondition(const std::string& field_name) const -> double;

    /**
     * @brief Gets the vector boundary condition value for a given field.
     * @param field_name The name of the field.
     * @return The vector value as Vector3d.
     * @throws std::runtime_error If the field is not found or is not a vector.
     */
    auto getVectorBoundaryCondition(const std::string& field_name) const -> Vector3d;

    /**
     * @brief Returns the const reference to face IDs associated with this patch.
     * @return const reference to vector of face IDs.
     */
    auto facesIds() const noexcept -> const std::vector<std::size_t>&;

    /**
     * @brief Returns the mutable reference to face IDs associated with this patch.
     * @return reference to vector of face IDs.
     */
    auto facesIds() noexcept -> std::vector<std::size_t>&;

    /**
     * @brief Checks if this patch is an empty patch.
     * @return true if the patch is empty, false otherwise.
     */
    auto isEmpty() const noexcept -> bool;

  private:
    auto getScalarBCSubfield(const std::string& name) const -> double;

    std::string _name;
    std::map<std::string, BoundaryCondition> _field_name_to_bc_map;
    std::vector<std::size_t> _faces_ids;
    bool _is_empty {false};
};

/**
 * @brief Represents the boundary configuration for a mesh, containing fields and boundary
 * patches.
 */
class MeshBoundary {
  public:
    /**
     * @brief Constructs a MeshBoundary by reading and parsing boundary condition files.
     * @param fields_path Path to the fields.json file.
     * @throws std::runtime_error If file cannot be opened, parsed, or contains invalid data.
     */
    explicit MeshBoundary(const std::filesystem::path& fields_path);

    /**
     * @brief Returns the vector of fields defined in the boundary file.
     * @return const reference to vector of FieldInfo.
     */
    auto fields() const noexcept -> const std::vector<FieldInfo>&;

    /**
     * @brief Returns the vector of boundary patches.
     * @return const reference to vector of BoundaryPatch.
     */
    auto patches() const noexcept -> const std::vector<BoundaryPatch>&;

  private:
    std::vector<FieldInfo> _fields;
    std::vector<BoundaryPatch> _boundary_patches;
};

} // namespace prism::mesh