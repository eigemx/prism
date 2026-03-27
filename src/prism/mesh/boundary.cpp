#include "boundary.h"

#include <fmt/base.h>
#include <fmt/format.h>

#include <fstream>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>

#include "prism/log.h"

using json = nlohmann::json;

namespace prism::mesh {

struct FieldBoundaryFile {
    std::string fieldName;
    std::map<std::string, BoundaryCondition> patchToBoundaryCondition;
};

/**
 * @brief Reads the entire contents of a file into a string.
 * @param path Path to the file to read.
 * @return The contents of the file as a string.
 */
auto fileToString(const std::filesystem::path& path) -> std::string {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error(
            fmt::format("prism::mesh::fileToString(): failed to open file {}", path.string()));
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    auto content = buffer.str();
    if (content.empty()) {
        throw std::runtime_error(fmt::format(
            "prism::mesh::fileToString(): boundary conditions file {} is empty!", path.string()));
    }
    return content;
}

/**
 * @brief Checks if a JSON string is valid.
 * @param json_data The JSON string to validate.
 * @return True if the JSON is valid, false otherwise.
 */
auto isGoodJson(const std::string& json_data) -> bool {
    return json::accept(json_data);
}

/**
 * @brief Checks if the JSON document contains a "fields" key.
 * @param doc The parsed JSON document.
 * @return True if "fields" key exists, false otherwise.
 */
auto containsFields(const json& doc) -> bool {
    return doc.contains(std::string("fields"));
}

/**
 * @brief Checks if the "fields" key in the JSON document is an array.
 * @param doc The parsed JSON document.
 * @return True if "fields" is an array, false otherwise.
 */
auto areFieldsAnArray(const json& doc) -> bool {
    return doc["fields"].is_array();
}

/**
 * @brief Checks if the JSON document contains a "patches" key.
 * @param doc The parsed JSON document.
 * @return True if "patches" key exists, false otherwise.
 */
auto containsPatches(const json& doc) -> bool {
    return doc.contains(std::string("patches"));
}

/**
 * @brief Checks if the "patches" key in the JSON document is an array.
 * @param doc The parsed JSON document.
 * @return True if "patches" is an array, false otherwise.
 */
auto arePatchesAnArray(const json& doc) -> bool {
    return doc["patches"].is_array();
}

/**
 * @brief Parses a single field entry from the JSON document.
 * @param field The JSON object representing a field.
 * @return A FieldInfo object containing the parsed field information.
 * @throws std::runtime_error If the field does not contain required "name" or "type" keys.
 */
auto parseField(const json& field) -> FieldInfo {
    std::optional<std::string> grad_scheme = std::nullopt;
    std::optional<std::vector<double>> units = std::nullopt;

    if (!field.contains("name") || !field.contains("type")) {
        throw std::runtime_error("prism::mesh::parseField(): Field must contain name and type");
    }

    auto name = field["name"].get<std::string>();
    auto type = field["type"].get<std::string>();

    if (field.contains("gradScheme")) {
        grad_scheme = field["gradScheme"].get<std::string>();
    }
    if (field.contains("units")) {
        units = field["units"].get<std::vector<double>>();
    }

    return {name, type, grad_scheme, units};
}

/**
 * @brief Parses all fields from the JSON document.
 * @param fields The JSON array of field objects.
 * @return A vector of FieldInfo objects.
 */
auto parseFields(const json& fields) -> std::vector<FieldInfo> {
    std::vector<FieldInfo> parsed_fields;
    for (const auto& field : fields) {
        parsed_fields.emplace_back(parseField(field));
    }
    return parsed_fields;
}

/**
 * @brief Checks if boundary condition files exist for all fields.
 * @param fields Vector of FieldInfo objects.
 * @param path Path to the fields.json file.
 * @return True if all field boundary files exist, false otherwise.
 */
auto fieldsFilesExist(const std::vector<FieldInfo>& fields, const std::filesystem::path& path)
    -> bool {
    for (const auto& field : fields) {
        auto file_name = fmt::format("{}.json", field.name());
        auto file_path = path.parent_path() / file_name;
        if (!std::filesystem::exists(file_path)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Reads and parses a single field's boundary condition file.
 * @param field_name The name of the field.
 * @param doc The parsed JSON document of the boundary file.
 * @return A FieldBoundaryFile containing the field name and patch-to-BC map.
 * @throws std::runtime_error If the file doesn't contain required "patches" array or patch
 * entries.
 */
auto readFieldsBoundaryFile(const std::string& field_name, const json& doc) -> FieldBoundaryFile {
    FieldBoundaryFile boundary_file;
    boundary_file.fieldName = field_name;

    if (!containsPatches(doc) || !arePatchesAnArray(doc)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): file `{}` either doesn't contain "
                        "`patches` field of `patches` entry is not an array",
                        field_name));
    }

    auto patches = doc["patches"];

    for (const auto& patch : patches) {
        if (!patch.contains("name") || !patch.contains("type")) {
            throw std::runtime_error(
                fmt::format("prism::mesh::readBoundaryFile(): file `{}` either doesn't contain "
                            "`name` or `type` field",
                            field_name));
        }
        auto name = patch["name"].get<std::string>();
        auto type = patch["type"].get<std::string>();

        // Check if patch has a value field
        if (patch.contains("value")) {
            auto value = patch["value"];
            if (value.is_number()) {
                auto bc = BoundaryCondition(
                    BoundaryConditionValueKind::Scalar, value.get<double>(), type);
                boundary_file.patchToBoundaryCondition.insert({name, bc});
            } else if (value.is_array()) {
                // make sure that the array is of size 3
                if (value.size() != 3) {
                    throw std::runtime_error(
                        fmt::format("prism::mesh::readBoundaryFile(): expected array of size 3 "
                                    "for field `{}` in patch `{}`",
                                    field_name,
                                    name));
                }
                auto value_array = value.get<std::vector<double>>();
                auto V = Vector3d {value_array[0], value_array[1], value_array[2]};
                auto bc = BoundaryCondition(BoundaryConditionValueKind::Vector, V, type);
                boundary_file.patchToBoundaryCondition.insert({name, bc});
            }
        } else {
            // patch does not have a value (empty, symmetry, outlet,...)
            auto bc = BoundaryCondition(BoundaryConditionValueKind::Nil, 0.0, type);
            boundary_file.patchToBoundaryCondition.insert({name, bc});
        }
    }
    return boundary_file;
}

/**
 * @brief Reads and parses all field boundary condition files.
 * @param path Path to the fields.json file.
 * @param fields Vector of FieldInfo objects.
 * @return A vector of FieldBoundaryFile objects.
 * @throws std::runtime_error If JSON parsing fails.
 */
auto readFieldsBoundaryFiles(const std::filesystem::path& path,
                             const std::vector<FieldInfo>& fields)
    -> std::vector<FieldBoundaryFile> {
    std::vector<FieldBoundaryFile> boundary_files;
    for (const auto& field : fields) {
        log::debug(
            "prism::mesh::readFieldsBoundaryFiles() : Reading boundary conditions file for field "
            "{}",
            field.name());
        auto file_name = fmt::format("{}.json", field.name());
        auto file_path = path.parent_path() / file_name;
        std::string json_data = fileToString(file_path);

        try {
            auto doc = json::parse(json_data);
            boundary_files.emplace_back(readFieldsBoundaryFile(field.name(), doc));
        } catch (json::parse_error& e) {
            throw std::runtime_error(fmt::format(
                "prism::mesh::readBoundaryFile(): Error parsing json file: {}", e.what()));
        }
    }

    return boundary_files;
}

/**
 * @brief Extracts unique boundary patches from all field boundary files.
 * @param boundary_files Vector of FieldBoundaryFile objects.
 * @return A vector of BoundaryPatch objects.
 */
auto extractBoundaryPatches(const std::vector<FieldBoundaryFile>& field_boundary_files)
    -> std::vector<BoundaryPatch> {
    auto patches = std::vector<BoundaryPatch>();

    for (const auto& [patch_name, _] : field_boundary_files.front().patchToBoundaryCondition) {
        std::map<std::string, BoundaryCondition> field_name_to_bc;
        for (const auto& boundary_file : field_boundary_files) {
            field_name_to_bc.insert(
                {boundary_file.fieldName, boundary_file.patchToBoundaryCondition.at(patch_name)});
        }

        patches.emplace_back(patch_name, field_name_to_bc);
    }
    return patches;
}

/**
 * @brief Parses boundary patches from field boundary files.
 * @param path Path to the fields.json file.
 * @param fields Vector of FieldInfo objects.
 * @return A vector of BoundaryPatch objects.
 * @throws std::runtime_error If field files don't exist or are inconsistent.
 */
auto parsePatches(const std::filesystem::path& path, const std::vector<FieldInfo>& fields)
    -> std::vector<BoundaryPatch> {
    // check if there is a json file for every field defined in fields.json file
    if (!fieldsFilesExist(fields, path)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Some field defined in `fields.json` "
                        "file are not present in a dedicated boundary conditions file"));
    }

    std::vector<FieldBoundaryFile> field_boundary_files = readFieldsBoundaryFiles(path, fields);

    if (field_boundary_files.empty()) {
        return {};
    }

    // make sure that boundary files are consistent
    auto patches_count = field_boundary_files.at(0).patchToBoundaryCondition.size();

    for (const auto& field_file : field_boundary_files) {
        if (field_file.patchToBoundaryCondition.size() != patches_count) {
            throw std::runtime_error(
                "prism::mesh::readBoundaryFile(): boundary conditions files are "
                "not consistent, please make sure that all boundary conditions "
                "files have the same number of patches");
        }
    }
    return extractBoundaryPatches(field_boundary_files);
}

/**
 * @brief Constructs a MeshBoundary by reading and parsing boundary condition files.
 * @param path Path to the fields.json file.
 * @throws std::runtime_error If file cannot be opened, parsed, or contains invalid data.
 */
MeshBoundary::MeshBoundary(const std::filesystem::path& path) {
    // read boundary file
    auto file = std::ifstream(path);

    if (!file) {
        throw std::runtime_error(fmt::format(
            "prism::mesh::readBoundaryFile(): Failed to open boundary conditions file `{}`",
            path.string()));
    }

    std::string json_data = fileToString(path);

    json doc;
    try {
        doc = json::parse(json_data);
    } catch (json::parse_error& e) {
        throw std::runtime_error(fmt::format(
            "prism::mesh::readBoundaryFile(): Error parsing json file: {}", e.what()));
    }

    if (!containsFields(doc) || !areFieldsAnArray(doc)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Couldn't find definition for fields in "
                        "boundary conditions file `{}`",
                        path.string()));
    }

    _fields = parseFields(doc["fields"]);
    _boundary_patches = parsePatches(path, _fields);
}

/**
 * @brief Returns the vector of fields defined in the boundary file.
 * @return const reference to vector of FieldInfo.
 */
auto MeshBoundary::fields() const noexcept -> const std::vector<FieldInfo>& {
    return _fields;
}

/**
 * @brief Returns the vector of boundary patches.
 * @return const reference to vector of BoundaryPatch.
 */
auto MeshBoundary::patches() const noexcept -> const std::vector<BoundaryPatch>& {
    return _boundary_patches;
}

/**
 * @brief Checks if a field name is a component field (e.g., U_x, U_y, U_z).
 * @param name The field name to check.
 * @return True if the field is a component of a vector field, false otherwise.
 */
auto isComponentField(const std::string& name) -> bool {
    if (name.size() < 2) {
        return false;
    }

    char last_char = name.back();
    if (last_char == 'x' || last_char == 'y' || last_char == 'z') {
        char second_last_char = *(name.rbegin() + 1);
        return second_last_char == '_';
    }

    return false;
}

/**
 * @brief Constructs a BoundaryCondition with the specified type, value, and kind string.
 * @param type The value kind (Nil, Scalar, or Vector).
 * @param value The boundary condition value.
 * @param bc_type_str The boundary condition type string (e.g., "fixedValue", "zeroGradient").
 */
BoundaryCondition::BoundaryCondition(BoundaryConditionValueKind type,
                                     BoundaryConditionValue value,
                                     std::string bc_type_str)
    : _value_kind(type), _value(std::move(value)), _kind_str(std::move(bc_type_str)) {}

/**
 * @brief Returns the value kind of this boundary condition.
 * @return The BoundaryConditionValueKind.
 */
auto BoundaryCondition::valueKind() const noexcept -> BoundaryConditionValueKind {
    return _value_kind;
}

/**
 * @brief Returns the value of this boundary condition.
 * @return const reference to the BoundaryConditionValue.
 */
auto BoundaryCondition::value() const noexcept -> const BoundaryConditionValue& {
    return _value;
}

/**
 * @brief Returns the kind string of this boundary condition.
 * @return const reference to the kind string.
 */
auto BoundaryCondition::kindString() const noexcept -> const std::string& {
    return _kind_str;
}

FieldInfo::FieldInfo(std::string name, std::string type)
    : _name(std::move(name)), _type(std::move(type)) {}

FieldInfo::FieldInfo(std::string name,
                     std::string type,
                     std::optional<std::string> grad_scheme,
                     std::optional<std::vector<double>> units)
    : _name(std::move(name)),
      _type(std::move(type)),
      _grad_scheme(std::move(grad_scheme)),
      _units(std::move(units)) {}

auto FieldInfo::name() const noexcept -> const std::string& {
    return _name;
}
auto FieldInfo::type() const noexcept -> const std::string& {
    return _type;
}
auto FieldInfo::gradScheme() const noexcept -> const std::optional<std::string>& {
    return _grad_scheme;
}
auto FieldInfo::units() const noexcept -> const std::optional<std::vector<double>>& {
    return _units;
}


/**
 * @brief Constructs a BoundaryPatch with the given name and field-to-BC map.
 * @param name The name of the boundary patch.
 * @param field_name_to_bc_map Map of field names to their boundary conditions.
 * @throws std::runtime_error If the field name to BC map is empty.
 */
BoundaryPatch::BoundaryPatch(std::string name,
                             std::map<std::string, BoundaryCondition> field_name_to_bc_map)
    : _name(std::move(name)), _field_name_to_bc_map(std::move(field_name_to_bc_map)) {
    if (_field_name_to_bc_map.empty()) {
        throw std::runtime_error(fmt::format(
            "BoundaryPatch constructor for patch {} was called with an empty field name to "
            "boundary condition map",
            _name));
    }

    const auto& it = _field_name_to_bc_map.begin();
    const auto& bc = it->second;

    if (bc.kindString() == "empty") {
        _is_empty = true;
    }
}

/**
 * @brief Returns the name of the boundary patch.
 * @return const reference to the patch name.
 */
auto BoundaryPatch::name() const noexcept -> const std::string& {
    return _name;
}

/**
 * @brief Checks if this patch is an empty patch.
 * @return true if the patch is empty, false otherwise.
 */
auto BoundaryPatch::isEmpty() const noexcept -> bool {
    return _is_empty;
}

/**
 * @brief Returns the const reference to face IDs associated with this patch.
 * @return const reference to vector of face IDs.
 */
auto BoundaryPatch::facesIds() const noexcept -> const std::vector<std::size_t>& {
    return _faces_ids;
}

/**
 * @brief Returns the mutable reference to face IDs associated with this patch.
 * @return reference to vector of face IDs.
 */
auto BoundaryPatch::facesIds() noexcept -> std::vector<std::size_t>& {
    return _faces_ids;
}

/**
 * @brief Gets the boundary condition for a given field name.
 * @param field_name The name of the field.
 * @return A const reference to the BoundaryCondition.
 * @throws std::runtime_error If the field is not found in the patch.
 */
auto BoundaryPatch::getBoundaryCondition(const std::string& field_name) const
    -> const BoundaryCondition& {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // search for the parent field, if exists
        if (isComponentField(field_name)) {
            it = _field_name_to_bc_map.find(field_name.substr(0, field_name.size() - 2));
        }
    }

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("prism::mesh::BoundaryPatch::getBoundaryCondition(): "
                        "Boundary patch `{}` does not have a field named `{}`",
                        _name,
                        field_name));
    }

    return it->second;
}

/**
 * @brief Gets the scalar boundary condition value for a given field.
 * @param field_name The name of the field.
 * @return The scalar value.
 * @throws std::runtime_error If the field is not found or is not a scalar.
 */
auto BoundaryPatch::getScalarBoundaryCondition(const std::string& field_name) const -> double {
    // in some cases we're dealing with a ScalarField that is a component of a parent
    // VectorField, such as when dealing with the x-component ScalarField of a velocity
    // VectorField U. in this case we won't find the boundary condition value for U_x as a
    // single scalar BC, but we can get it from the first component of the vector BC value of
    // U.
    if (isComponentField(field_name)) {
        return getScalarBCSubfield(field_name);
    }

    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("prism::mesh::BoundaryPatch::getScalarBoundaryCondition(): "
                        "Boundary patch '{}' does not have a field named '{}'",
                        _name,
                        field_name));
    }

    if (it->second.valueKind() != BoundaryConditionValueKind::Scalar) {
        // field is not a scalar
        throw std::runtime_error(
            fmt::format("prism::mesh::BoundaryPatch::getScalarBoundaryCondition(): "
                        "Boundary patch '{}' field '{}' is not a scalar",
                        _name,
                        field_name));
    }

    return std::get<double>(it->second.value());
}


/**
 * @brief Gets the vector boundary condition value for a given field.
 * @param field_name The name of the field.
 * @return The vector value as Vector3d.
 * @throws std::runtime_error If the field is not found or is not a vector.
 */
auto BoundaryPatch::getVectorBoundaryCondition(const std::string& field_name) const -> Vector3d {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("prism::mesh::BoundaryPatch::getVectorBoundaryCondition(): "
                        "Boundary patch '{}' does not have a field named '{}'",
                        _name,
                        field_name));
    }

    if (it->second.valueKind() != BoundaryConditionValueKind::Vector) {
        // field is not a vector
        throw std::runtime_error(
            fmt::format("prism::mesh::BoundaryPatch::getVectorBoundaryCondition(): "
                        "Boundary patch '{}' field '{}' is not a vector",
                        _name,
                        field_name));
    }

    return std::get<Vector3d>(it->second.value());
}

/**
 * @brief Gets a scalar component from a vector boundary condition subfield.
 * @param name The component field name (e.g., "U_x", "U_y", "U_z").
 * @return The scalar component value.
 * @throws std::runtime_error If the field name doesn't end with a valid cartesian component.
 */
auto BoundaryPatch::getScalarBCSubfield(const std::string& name) const -> double {
    const auto& vec_value = getVectorBoundaryCondition(name.substr(0, name.size() - 2));

    switch (name.back()) {
        case 'x': return vec_value[0];
        case 'y': return vec_value[1];
        case 'z': return vec_value[2];
        default: {
            throw std::runtime_error(
                "prism::mesh::BoundaryPatch::getScalarBCSubfield() was given a field with a "
                "name that does not end with a valid cartesian component (aka x, y or z).");
        }
    }
}


} // namespace prism::mesh
