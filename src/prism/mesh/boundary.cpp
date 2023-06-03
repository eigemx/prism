// TODO: boundary file reading results in segmenation fault in case boundary file is not as per the
// expected format.
#include "boundary.h"

#include <toml++/toml.h>

#include <fstream>
#include <iostream>
#include <string_view>
#include <unordered_map>

#include "../print.h"

namespace prism::mesh {
auto boundary_type_str_to_enum(std::string_view type) -> BoundaryPatchType;
auto parse_boundary_patch(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch;
auto parse_nested_boundary_conditions(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch;
auto parse_field_boundary_condition(const toml::table& table,
                                    const std::string& patch_name,
                                    const std::string& field_name) -> BoundaryCondition;


auto boundary_type_str_to_enum(std::string_view type) -> BoundaryPatchType {
    const auto static bc_type_map = std::unordered_map<std::string_view, BoundaryPatchType> {
        {"fixed", BoundaryPatchType::Fixed},
        {"inlet", BoundaryPatchType::Inlet},
        {"outlet", BoundaryPatchType::Outlet},
        {"gradient", BoundaryPatchType::FixedGradient},
        {"symmetry", BoundaryPatchType::Symmetry},
        {"empty", BoundaryPatchType::Empty},
    };

    auto it = bc_type_map.find(type);

    if (it == bc_type_map.end()) {
        return BoundaryPatchType::Unknown;
    }

    return it->second;
}

auto parse_boundary_patch(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch {
    return parse_nested_boundary_conditions(table, patch_name);
}

auto parse_nested_boundary_conditions(const toml::table& table, const std::string& patch_name)
    -> BoundaryPatch {
    // for each subtable, get its type and data
    std::map<std::string, BoundaryCondition> field_name_to_bc_map;
    const auto& patch_table = *(table[patch_name].as_table());

    for (const auto& [field_name_key, field_table] : patch_table) {
        const auto& field_name = std::string(field_name_key.str());
        auto field_bc = parse_field_boundary_condition(table, patch_name, field_name);

        field_name_to_bc_map.insert({field_name, field_bc});

        if (field_bc.type() == BoundaryConditionValueType::Vector) {
            // to make it easier to access the x, y, z boundary conditions for a vector field
            auto vec_value = std::get<Vector3d>(field_bc.data());
            field_name_to_bc_map.insert({field_name + "_x",
                                         BoundaryCondition {
                                             BoundaryConditionValueType::Scalar,
                                             vec_value.x(),
                                             field_bc.patch_type(),
                                         }});
            field_name_to_bc_map.insert({field_name + "_y",
                                         BoundaryCondition {
                                             BoundaryConditionValueType::Scalar,
                                             vec_value.y(),
                                             field_bc.patch_type(),
                                         }});
            field_name_to_bc_map.insert({field_name + "_z",
                                         BoundaryCondition {
                                             BoundaryConditionValueType::Scalar,
                                             vec_value.z(),
                                             field_bc.patch_type(),
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
            format("boundary.cpp parse_field_boundary_condition(): "
                   "Boundary patch '{}' field '{}' does not have a type or value.",
                   patch_name,
                   field_name));
    }

    auto bc_type_str = field_table["type"].value<std::string_view>().value();
    auto bc_type = boundary_type_str_to_enum(bc_type_str);


    if (!field_table.contains("value")) {
        return BoundaryCondition {
            BoundaryConditionValueType::Nil, BoundaryConditionData {}, bc_type};
    }

    auto bc_value = field_table["value"];

    if (bc_value.is_number() || bc_value.is_floating_point()) {
        double value = bc_value.value<double>().value();
        return BoundaryCondition {BoundaryConditionValueType::Scalar, value, bc_type};
    }

    if (bc_value.is_array()) {
        const auto& array = bc_value.as_array();

        if (array->size() != 3) {
            throw std::runtime_error(
                format("boundary.cpp parse_field_boundary_condition(): "
                       "Array value for field '{}' for patch '{}' is not a 3D vector",
                       field_name,
                       patch_name));
        }

        return {BoundaryConditionValueType::Vector,
                Vector3d {
                    array->at(0).value<double>().value(),
                    array->at(1).value<double>().value(),
                    array->at(2).value<double>().value(),
                },
                bc_type};
    }

    throw std::runtime_error(
        format("boundary.cpp parse_field_boundary_condition(): "
               "Boundary patch '{}' field '{}' has an invalid type or value.",
               patch_name,
               field_name));
}


auto read_boundary_data_file(const std::filesystem::path& path,
                             const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch> {
    std::vector<BoundaryPatch> boundary_patches;
    boundary_patches.reserve(boundary_names.size());

    auto fstream {std::ifstream {path}};

    // throw if file doesn't exist
    if (!fstream) {
        throw std::runtime_error(
            format("Failed to open boundary conditions file '{}'", path.string()));
    }

    toml::table doc;

    // parse boundary TOML file
    try {
        doc = toml::parse(fstream);
    } catch (const toml::parse_error& e) {
        // file cannot be parsed
        throw std::runtime_error(
            format("Failed to parse boundary condition file: '{}', complete error: {}",
                   path.string(),
                   e.what()));
    }

    // for each defined boundary, get its relevant BoundaryData object
    for (const auto& bname : boundary_names) {
        auto table {doc[bname.data()]};

        if (!table) {
            throw std::runtime_error(
                format("Couldn't find definition for boundary patch '{}' in boundary "
                       "conditions file '{}'",
                       bname,
                       path.string()));
        }

        boundary_patches.emplace_back(parse_boundary_patch(doc, std::string(bname)));
    }

    return boundary_patches;
}

auto BoundaryPatch::get_scalar_bc(const std::string& field_name) const -> double {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            format("BoundaryPatch::get_scalar_bc(): "
                   "Boundary patch '{}' does not have a field named '{}'",
                   _name,
                   field_name));
    }

    if (it->second.type() != BoundaryConditionValueType::Scalar) {
        // field is not a scalar
        throw std::runtime_error(
            format("BoundaryPatch::get_scalar_bc(): "
                   "Boundary patch '{}' field '{}' is not a scalar",
                   _name,
                   field_name));
    }

    return std::get<double>(it->second.data());
}

auto BoundaryPatch::get_vector_bc(const std::string& field_name) const -> Vector3d {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            format("BoundaryPatch::get_vector_bc(): "
                   "Boundary patch '{}' does not have a field named '{}'",
                   _name,
                   field_name));
    }

    if (it->second.type() != BoundaryConditionValueType::Vector) {
        // field is not a vector
        throw std::runtime_error(
            format("BoundaryPatch::get_vector_bc(): "
                   "Boundary patch '{}' field '{}' is not a vector",
                   _name,
                   field_name));
    }

    return std::get<Vector3d>(it->second.data());
}

auto BoundaryPatch::get_bc(const std::string& field_name) const -> const BoundaryCondition& {
    // Search for the field name in the boundary patch
    auto it = _field_name_to_bc_map.find(field_name);

    if (it == _field_name_to_bc_map.end()) {
        // field not found
        throw std::runtime_error(
            format("BoundaryPatch::get_bc(): "
                   "Boundary patch '{}' does not have a field named '{}'",
                   _name,
                   field_name));
    }

    return it->second;
}

} // namespace prism::mesh