#include "boundary.h"

#include <toml++/toml.h>

#include <fstream>
#include <iostream>
#include <string_view>
#include <unordered_map>

#include "../print.h"

namespace prism::mesh {

auto BoundaryPatch::get_scalar_bc(const std::string& name) const -> double {
    for (const auto& boundary_condition : _bcs) {
        // In case we are looking for a scalar boundary condition, like temperature.
        if (boundary_condition.name() == name &&
            boundary_condition.type() == BoundaryConditionType::Scalar) {
            return std::get<double>(boundary_condition.data());
        }

        // In case we are looking for a scalar boundary condition, inside a vector boundary condition, like velocity.
        // For example, we are looking for x-component of inlet velocity.
        // in that case `name` parameter will end with `_x` and we need to remove that `_x` suffix.
        // and then check if the name matches and also the boundary condition type is vector.
        // if yes, then we return the x or y or z component of the vector (depending on the suffix).
        if (boundary_condition.name() == name.substr(0, name.size() - 2) &&
            boundary_condition.type() == BoundaryConditionType::Vector) {
            auto vector = std::get<Vector3d>(boundary_condition.data());
            auto suffix = std::string_view(name).substr(name.size() - 2, 2);

            if (suffix == "_x") {
                return vector.x();
            }

            if (suffix == "_y") {
                return vector.y();
            }

            if (suffix == "_z") {
                return vector.z();
            }
        }
    }

    throw std::runtime_error(
        format("mesh::BoundaryPatch::get_scalar_bc(): "
               "Boundary patch '{}' does not have boundary condition '{}'",
               _name,
               name));
}

auto BoundaryPatch::get_vector_bc(const std::string& name) const -> Vector3d {
    for (const auto& attr : _bcs) {
        if (attr.name() == name && attr.type() == BoundaryConditionType::Vector) {
            return std::get<Vector3d>(attr.data());
        }
    }

    throw std::runtime_error(
        format("mesh::BoundaryPatch::get_vector_bc(): "
               "Boundary patch '{}' does not have boundary condition '{}'",
               _name,
               name));
}

auto boundary_type_str_to_enum(std::string_view type) -> BoundaryPatchType {
    const auto static bc_type_map = std::unordered_map<std::string_view, BoundaryPatchType> {
        {"wall", BoundaryPatchType::Fixed},
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

auto parse_boundary_conditions(const toml::table& table, std::string_view bname)
    -> BoundaryConditions {
    BoundaryConditions boundary_conditions;

    // read all toml nodes in table, if of type double or array and of size 3, add to boundary_conditions
    const auto& toml_nodes = table[bname].as_table();
    for (auto&& [key, value] : *toml_nodes) {
        if (key == "type") {
            continue;
        }

        if (value.is_number() || value.is_floating_point()) {
            double attr_float_val = value.value<double>().value();
            boundary_conditions.emplace_back(
                std::string(key), BoundaryConditionType::Scalar, attr_float_val);
        }

        else if (value.is_array()) {
            const auto& array = value.as_array();

            if (array->size() != 3) {
                throw std::runtime_error(
                    format("Boundary condition '{}' for patch '{}' is not a 3D vector",
                           key.str(),
                           bname));
            }

            boundary_conditions.emplace_back(std::string(key),
                                             BoundaryConditionType::Vector,
                                             Vector3d {
                                                 array->at(0).value<double>().value(),
                                                 array->at(1).value<double>().value(),
                                                 array->at(2).value<double>().value(),
                                             });
        }

        else {
            throw std::runtime_error(
                format("Boundary condition '{}' for patch '{}' is not a scalar or vector",
                       key.str(),
                       bname));
        }
    }

    return boundary_conditions;
}

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> std::vector<BoundaryPatch> {
    std::vector<BoundaryPatch> bcs;
    bcs.reserve(boundary_names.size());

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

        auto type {doc[bname]["type"].value<std::string_view>()};

        if (!type) {
            throw std::runtime_error(format(
                "Boundary conditions for patch '{}' does not have a type in boundary condition "
                "file '{}'",
                bname,
                path.string()));
        }

        auto boundary_patch_type = boundary_type_str_to_enum(type.value());
        auto bp_attributes = parse_boundary_conditions(doc, bname);

        bcs.emplace_back(std::string(bname), std::move(bp_attributes), boundary_patch_type);
    }

    return bcs;
}
} // namespace prism::mesh