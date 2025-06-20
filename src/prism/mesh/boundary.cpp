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

auto fileToString(const std::filesystem::path& path) -> std::string {
    std::ifstream file(path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

auto isGoodJson(const std::string& json_data) -> bool {
    return json::accept(json_data);
}

auto containsFields(const json& doc) -> bool {
    return doc.contains(std::string("fields"));
}

auto areFieldsAnArray(const json& doc) -> bool {
    return doc["fields"].is_array();
}

auto containsPatches(const json& doc) -> bool {
    return doc.contains(std::string("patches"));
}

auto arePatchesAnArray(const json& doc) -> bool {
    return doc["patches"].is_array();
}

auto parseField(const json& field) -> FieldInfo {
    std::optional<std::string> grad_scheme = std::nullopt;
    std::optional<std::vector<double>> units = std::nullopt;

    if (!field.contains("name") || !field.contains("type")) {
        throw std::runtime_error("Field must contain name and type");
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

auto parseFields(const json& fields) -> std::vector<FieldInfo> {
    std::vector<FieldInfo> parsed_fields;
    for (const auto& field : fields) {
        parsed_fields.emplace_back(parseField(field));
    }
    return parsed_fields;
}

auto fieldsFilesExist(const std::vector<FieldInfo>& fields,
                      const std::filesystem::path& path) -> bool {
    for (const auto& field : fields) {
        auto file_name = fmt::format("{}.json", field.name());
        auto file_path = path.parent_path() / file_name;
        if (!std::filesystem::exists(file_path)) {
            return false;
        }
    }
    return true;
}

auto readFieldsBoundaryFile(const std::string& field_name, const json& doc) -> FieldBoundaryFile {
    FieldBoundaryFile boundary_file;
    boundary_file.fieldName = field_name;

    if (!containsPatches(doc) || !arePatchesAnArray(doc)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): file `{}` either don't contain "
                        "`patches` field of `patches` is not an array",
                        field_name));
    }

    auto patches = doc["patches"];

    for (const auto& patch : patches) {
        if (!patch.contains("name") || !patch.contains("type")) {
            throw std::runtime_error(
                fmt::format("prism::mesh::readBoundaryFile(): file `{}` either don't contain "
                            "`name` or `type` field",
                            field_name));
        }
        auto name = patch["name"].get<std::string>();
        auto type = patch["type"].get<std::string>();

        if (patch.contains("value")) {
            auto value = patch["value"];
            if (value.is_number()) {
                auto bc = BoundaryCondition(
                    BoundaryConditionValueKind::Scalar, value.get<double>(), type);
                boundary_file.patchToBoundaryCondition.insert({name, bc});
            }

            if (value.is_array()) {
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
        }
        // patch do not have a value (empty, symmetry, outlet,...)
        auto bc = BoundaryCondition(BoundaryConditionValueKind::Nil, 0.0, type);
        boundary_file.patchToBoundaryCondition.insert({name, bc});
    }
    return boundary_file;
}

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

        if (!isGoodJson(json_data)) {
            try {
                json Doc {json::parse(json_data)};
            } catch (json::parse_error& e) {
                throw std::runtime_error(fmt::format(
                    "prism::mesh::readBoundaryFile(): Error parsing json file: {}", e.what()));
            }
        }

        auto doc = json::parse(json_data);
        boundary_files.emplace_back(readFieldsBoundaryFile(field.name(), doc));
    }

    return boundary_files;
}

auto extractBoundaryPatches(const std::vector<FieldBoundaryFile>& boundary_files)
    -> std::vector<BoundaryPatch> {
    // each boundary file contains a map of patch name to boundary condition, we need to extract
    // the patch names (keys  of the map) to a std::vector
    std::set<std::string> patch_names;
    for (const auto& boundary_file : boundary_files) {
        for (const auto& [name, _] : boundary_file.patchToBoundaryCondition) {
            patch_names.insert(name);
        }
    }

    auto patches = std::vector<BoundaryPatch>();
    for (const auto& patch_name : patch_names) {
        std::map<std::string, BoundaryCondition> field_to_bc;
        for (const auto& boundary_file : boundary_files) {
            field_to_bc.insert(
                {boundary_file.fieldName, boundary_file.patchToBoundaryCondition.at(patch_name)});
        }

        patches.emplace_back(patch_name, field_to_bc);
    }
    return patches;
}

auto parsePatches(const std::filesystem::path& path,
                  const std::vector<FieldInfo>& fields) -> std::vector<BoundaryPatch> {
    // check if there is a json file for every field
    if (!fieldsFilesExist(fields, path)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Some field defined in `fields.json` "
                        "file are not present in a dedicated boundary conditions file"));
    }

    auto boundary_files = readFieldsBoundaryFiles(path, fields);

    // make sure that boundary files are consistent
    auto patches_count = boundary_files.at(0).patchToBoundaryCondition.size();

    for (const auto& boundary_file : boundary_files) {
        if (boundary_file.patchToBoundaryCondition.size() != patches_count) {
            throw std::runtime_error(
                "prism::mesh::readBoundaryFile(): boundary conditions files are "
                "not consistent, please make sure that all boundary conditions "
                "files have the same number of patches");
        }
    }
    return extractBoundaryPatches(boundary_files);
}

MeshBoundary::MeshBoundary(const std::filesystem::path& path) {
    // read boundary file
    auto file = std::ifstream(path);

    if (!file) {
        throw std::runtime_error(fmt::format(
            "prism::mesh::readBoundaryFile(): Failed to open boundary conditions file `{}`",
            path.string()));
    }

    std::string json_data = fileToString(path);

    if (!isGoodJson(json_data)) {
        try {
            json Doc {json::parse(json_data)};
        } catch (json::parse_error& e) {
            throw std::runtime_error(fmt::format(
                "prism::mesh::readBoundaryFile(): Error parsing json file: {}", e.what()));
        }
    }

    auto doc = json::parse(json_data);

    if (!containsFields(doc) || !areFieldsAnArray(doc)) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Couldn't find definition for fields in "
                        "boundary conditions file `{}`",
                        path.string()));
    }

    _fields = parseFields(doc["fields"]);
    _boundary_patches = parsePatches(path, _fields);
}

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

BoundaryCondition::BoundaryCondition(BoundaryConditionValueKind type,
                                     BoundaryConditionValue value,
                                     std::string bc_type_str)
    : _value_kind(type), _value(std::move(value)), _kind_str(std::move(bc_type_str)) {}


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