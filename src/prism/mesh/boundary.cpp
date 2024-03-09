// TODO: boundary file reading results in segmenation fault in case boundary file is not as per
// the expected format.
#include "boundary.h"

#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <toml++/toml.h>

#include <stdexcept>
#include <string_view>
#include <unordered_map>

// TODO: We're using name.substr(0, name.size() - 2) many times to get the parent field name
// it's better to wrap this in a little inline function

namespace prism::mesh {
// Parsing boundary file functions
auto inline boundary_type_str_to_enum(std::string_view type) -> BoundaryConditionKind;
auto inline parse_boundary_patch(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch;
auto inline parse_nested_boundary_conditions(const toml::table& table,
                                             const std::string& patch_name) -> BoundaryPatch;
auto inline parse_field_boundary_condition(const toml::table& table,
                                           const std::string& patch_name,
                                           const std::string& field_name) -> BoundaryCondition;
auto is_component_field(const std::string& name) -> bool;

auto boundary_type_str_to_enum(std::string_view type) -> BoundaryConditionKind {
    const auto static bc_type_map = std::unordered_map<std::string_view, BoundaryConditionKind> {
        {"fixed", BoundaryConditionKind::Fixed},
        // TODO: Should this ve velocityInlet?
        {"inlet", BoundaryConditionKind::VelocityInlet},
        {"outlet", BoundaryConditionKind::Outlet},
        {"gradient", BoundaryConditionKind::FixedGradient},
        {"symmetry", BoundaryConditionKind::Symmetry},
        {"empty", BoundaryConditionKind::Empty},
        // TODO: Implement remaining boundary condition kinds
    };

    auto it = bc_type_map.find(type);

    if (it == bc_type_map.end()) {
        spdlog::warn("Boundary conditions file contains an unknown boundary patch type {}", type);
        return BoundaryConditionKind::Unknown;
    }

    return it->second;
}

auto parse_boundary_patch(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch {
    return parse_nested_boundary_conditions(table, patch_name);
}

auto parse_nested_boundary_conditions(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch {
    // for each sub-table, get its type and data
    std::map<std::string, BoundaryCondition> field_name_to_bc_map;
    const auto& patch_table = *(table[patch_name].as_table());

    for (const auto& [field_name_key, field_table] : patch_table) {
        const auto& field_name = std::string(field_name_key.str());
        auto field_bc = parse_field_boundary_condition(table, patch_name, field_name);

        field_name_to_bc_map.insert({field_name, field_bc});

        if (field_bc.value_kind() == BoundaryConditionValueKind::Vector) {
            // to make it easier to access the x, y, z boundary conditions for a vector field
            auto vec_value = std::get<Vector3d>(field_bc.value());
            field_name_to_bc_map.insert({field_name + "_x",
                                         BoundaryCondition {
                                             BoundaryConditionValueKind::Scalar,
                                             vec_value.x(),
                                             field_bc.kind(),
                                             field_bc.kind_string(),
                                         }});
            field_name_to_bc_map.insert({field_name + "_y",
                                         BoundaryCondition {
                                             BoundaryConditionValueKind::Scalar,
                                             vec_value.y(),
                                             field_bc.kind(),
                                             field_bc.kind_string(),
                                         }});
            field_name_to_bc_map.insert({field_name + "_z",
                                         BoundaryCondition {
                                             BoundaryConditionValueKind::Scalar,
                                             vec_value.z(),
                                             field_bc.kind(),
                                             field_bc.kind_string(),
                                         }});
        }
    }

    return BoundaryPatch {patch_name, field_name_to_bc_map};
}

auto parse_field_boundary_condition(const toml::table& table,
                                    const std::string& patch_name,
                                    const std::string& field_name) -> BoundaryCondition {
    const auto& field_table = *(table[patch_name][field_name].as_table());

    if (!field_table.contains("type")) {
        throw std::runtime_error(
            fmt::format("mesh/boundary.cpp::parse_field_boundary_condition(): "
                        "Boundary patch '{}' field '{}' does not have a type or value.",
                        patch_name,
                        field_name));
    }

    auto bc_type_str = field_table["type"].value<std::string_view>().value();
    auto bc_type = boundary_type_str_to_enum(bc_type_str);


    if (!field_table.contains("value")) {
        return BoundaryCondition {BoundaryConditionValueKind::Nil,
                                  BoundaryConditionValue {},
                                  bc_type,
                                  std::string(bc_type_str)};
    }

    auto bc_value = field_table["value"];

    if (bc_value.is_number() || bc_value.is_floating_point()) {
        double value = bc_value.value<double>().value();
        return BoundaryCondition {
            BoundaryConditionValueKind::Scalar, value, bc_type, std::string(bc_type_str)};
    }

    if (bc_value.is_array()) {
        const auto& array = bc_value.as_array();

        if (array->size() != 3) {
            throw std::runtime_error(
                fmt::format("mesh/boundary.cpp::parse_field_boundary_condition(): "
                            "Array value for field '{}' for patch '{}' is not a 3D vector",
                            field_name,
                            patch_name));
        }

        return {BoundaryConditionValueKind::Vector,
                Vector3d {
                    array->at(0).value<double>().value(),
                    array->at(1).value<double>().value(),
                    array->at(2).value<double>().value(),
                },
                bc_type,
                std::string(bc_type_str)};
    }

    throw std::runtime_error(
        fmt::format("mesh/boundary.cpp::parse_field_boundary_condition(): "
                    "Boundary patch '{}' field '{}' has an invalid type or value.",
                    patch_name,
                    field_name));
}


auto read_boundary_file(const std::filesystem::path& path,
                        const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch> {
    std::vector<BoundaryPatch> boundary_patches;
    boundary_patches.reserve(boundary_names.size());

    auto fstream {std::ifstream {path}};

    // throw if file doesn't exist
    if (!fstream) {
        throw std::runtime_error(
            fmt::format("Failed to open boundary conditions file '{}'", path.string()));
    }

    toml::table doc;

    // parse boundary TOML file
    try {
        doc = toml::parse(fstream);
    } catch (const toml::parse_error& e) {
        // file cannot be parsed
        throw std::runtime_error(
            fmt::format("Failed to parse boundary condition file: '{}', complete error: {}",
                        path.string(),
                        e.what()));
    }

    // for each defined boundary, get its relevant BoundaryData object
    for (const auto& bname : boundary_names) {
        auto table {doc[bname.data()]};

        if (!table) {
            throw std::runtime_error(
                fmt::format("Couldn't find definition for boundary patch '{}' in boundary "
                            "conditions file '{}'",
                            bname,
                            path.string()));
        }

        boundary_patches.emplace_back(parse_boundary_patch(doc, std::string(bname)));
    }

    return boundary_patches;
}


auto is_component_field(const std::string& name) -> bool {
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


auto BoundaryPatch::get_bc(const std::string& field_name) const -> const BoundaryCondition& {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // search for the parent field, if exists
        if (is_component_field(field_name)) {
            it = _field_name_to_bc_map.find(field_name.substr(0, field_name.size() - 2));
        }
    }

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("mesh/boundary.cpp::BoundaryPatch::get_bc(): "
                        "Boundary patch '{}' does not have a field named '{}'",
                        _name,
                        field_name));
    }

    return it->second;
}


auto BoundaryPatch::get_scalar_bc(const std::string& field_name) const -> double {
    // in some cases we're dealing with a ScalarField that is a component of a parent VectorField,
    // such as when dealing with the x-component ScalarField of a velocity VectorField U.
    // in this case we won't find the boundary condition value for U_x as a single scalar BC,
    // but we can get it from the first component of the vector BC value of U.
    if (is_component_field(field_name)) {
        return get_scalar_bc_subfield(field_name);
    }

    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("BoundaryPatch::get_scalar_bc(): "
                        "Boundary patch '{}' does not have a field named '{}'",
                        _name,
                        field_name));
    }

    if (it->second.value_kind() != BoundaryConditionValueKind::Scalar) {
        // field is not a scalar
        throw std::runtime_error(
            fmt::format("BoundaryPatch::get_scalar_bc(): "
                        "Boundary patch '{}' field '{}' is not a scalar",
                        _name,
                        field_name));
    }

    return std::get<double>(it->second.value());
}


auto BoundaryPatch::get_vector_bc(const std::string& field_name) const -> Vector3d {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            fmt::format("BoundaryPatch::get_vector_bc(): "
                        "Boundary patch '{}' does not have a field named '{}'",
                        _name,
                        field_name));
    }

    if (it->second.value_kind() != BoundaryConditionValueKind::Vector) {
        // field is not a vector
        throw std::runtime_error(
            fmt::format("BoundaryPatch::get_vector_bc(): "
                        "Boundary patch '{}' field '{}' is not a vector",
                        _name,
                        field_name));
    }

    return std::get<Vector3d>(it->second.value());
}

auto BoundaryPatch::get_scalar_bc_subfield(const std::string& name) const -> double {
    const auto& vec_value = get_vector_bc(name.substr(0, name.size() - 2));

    switch (name.back()) {
        case 'x':
            return vec_value[0];
        case 'y':
            return vec_value[1];
        case 'z':
            return vec_value[2];
        default: {
            throw std::runtime_error(
                "BoundaryPatch::get_scalar_bc_subfield was given a field with a name that does "
                "not in with a valid Cartesian componenet.");
        }
    }

    // It should be impossible to reach this
    return 0.0;
}

} // namespace prism::mesh