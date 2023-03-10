#include "boundary.h"

#include <fmt/format.h>
#include <toml++/toml.h>

#include <fstream>
#include <iostream>
#include <string_view>
#include <unordered_map>

namespace prism::mesh {

auto boundary_type_str_to_enum(std::string_view type) -> BoundaryPatchType {
    const auto static bc_type_map = std::unordered_map<std::string_view, BoundaryPatchType> {
        {"wall", BoundaryPatchType::Wall},
        {"inlet", BoundaryPatchType::Inlet},
        {"outlet", BoundaryPatchType::Outlet},
        {"gradient", BoundaryPatchType::Gradient},
        {"symmetry", BoundaryPatchType::Symmetry},
        {"empty", BoundaryPatchType::Empty},
    };

    auto it = bc_type_map.find(type);

    if (it == bc_type_map.end()) {
        return BoundaryPatchType::Unknown;
    }

    return it->second;
}

template <typename T = const toml::node_view<const toml::node>>
inline void check_velocity_or_temp(T velocity,
                                   T temperature,
                                   std::string_view bname,
                                   std::string_view btype) {
    if (!velocity && !temperature) {
        throw std::runtime_error(
            fmt::format("Boundary '{}' is defined as a '{}' and has no velocity or temperature",
                        bname,
                        btype));
    }
}

inline auto get_real_value(const toml::node_view<const toml::node>& node,
                           std::string_view bname,
                           std::string_view btype,
                           std::string_view value_type) -> double {
    if (!node.is_floating_point() && !node.is_number()) {
        throw std::runtime_error(fmt::format(
            "{} boundary '{}' has a {} that is not a numeric value", btype, bname, value_type));
    }

    return node.value<double>().value();
}

inline auto get_velocity(const toml::node_view<const toml::node>& velocity,
                         std::string_view bname,
                         std::string_view btype) -> Vector3d {
    // check if velocity is a vector
    if (!velocity.is_array() || velocity.as_array()->size() != 3) {
        throw std::runtime_error(
            fmt::format("{} boundary '{}' has a velocity that is not a 3D vector", btype, bname));
    }

    // check if velocity is a vector of numbers
    for (std::size_t i = 0; i < 3; ++i) {
        if (!velocity[i].is_floating_point() && !velocity[i].is_number()) {
            throw std::runtime_error(
                fmt::format("{} boundary '{}' has a velocity that is not a 3D vector of numbers",
                            btype,
                            bname));
        }
    }

    return Vector3d {velocity[0].value<double>().value(),
                     velocity[1].value<double>().value(),
                     velocity[2].value<double>().value()};
}


auto parse_wall(const toml::table& doc, std::string_view bname) -> WallBoundaryData {
    auto velocity = doc[bname]["velocity"];
    auto temperature = doc[bname]["temperature"];

    check_velocity_or_temp(velocity, temperature, bname, "Wall");

    std::optional<Vector3d> velocity_vec {std::nullopt};
    std::optional<double> temperature_val {std::nullopt};

    if (velocity) {
        velocity_vec = get_velocity(velocity, bname, "Wall");
    }

    if (temperature) {
        temperature_val = get_real_value(temperature, bname, "Wall", "temperature");
    }

    return WallBoundaryData {velocity_vec, temperature_val};
}

auto parse_inlet(const toml::table& doc, std::string_view bname) -> InletBoundaryData {
    auto velocity = doc[bname]["velocity"];
    auto temperature = doc[bname]["temperature"];

    check_velocity_or_temp(velocity, temperature, bname, "Inlet");

    std::optional<Vector3d> velocity_vec {std::nullopt};
    std::optional<double> temperature_val {std::nullopt};

    if (velocity) {
        velocity_vec = get_velocity(velocity, bname, "Inlet");
    }

    if (temperature) {
        temperature_val = get_real_value(temperature, bname, "Inlet", "temperature");
    }

    return InletBoundaryData {velocity_vec, temperature_val};
}

auto parse_outlet(const toml::table& doc, std::string_view bname) -> OutletBoundaryData {
    auto pressure = doc[bname]["pressure"];

    if (!pressure) {
        throw std::runtime_error(fmt::format(
            "Boundary '{}' is defined as an outlet but has no pressure value", bname));
    }

    double pressure_val = get_real_value(pressure, bname, "Outlet", "pressure");

    return OutletBoundaryData {pressure_val};
}

auto parse_boundary_data(const toml::table& table,
                         std::string_view bname,
                         BoundaryPatchType bc_type) -> BoundaryData {
    switch (bc_type) {
        case BoundaryPatchType::Wall:
            return parse_wall(table, bname);

        case BoundaryPatchType::Inlet:
            return parse_inlet(table, bname);

        case BoundaryPatchType::Outlet:
            return parse_outlet(table, bname);

        case BoundaryPatchType::Gradient:
            throw std::runtime_error("Gradient boundary condition is not implemented");

        case BoundaryPatchType::Symmetry:
            return SymmetryBoundaryData {};

        case BoundaryPatchType::Empty:
            return EmptyBoundaryData {};

        case BoundaryPatchType::Unknown:
            throw std::runtime_error("Unknown boundary condition type");
    }
}

auto BoundaryPatch::infer_boundary_type(const BoundaryData& data) -> BoundaryPatchType {
    if (std::holds_alternative<WallBoundaryData>(data)) {
        return BoundaryPatchType::Wall;
    }
    if (std::holds_alternative<InletBoundaryData>(data)) {
        return BoundaryPatchType::Inlet;
    }
    if (std::holds_alternative<OutletBoundaryData>(data)) {
        return BoundaryPatchType::Outlet;
    }
    if (std::holds_alternative<GradientBoundaryData>(data)) {
        return BoundaryPatchType::Gradient;
    }
    if (std::holds_alternative<SymmetryBoundaryData>(data)) {
        return BoundaryPatchType::Symmetry;
    }
    if (std::holds_alternative<EmptyBoundaryData>(data)) {
        return BoundaryPatchType::Empty;
    }
    throw std::runtime_error("Non-implemented boundary type");
}

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> BoundaryPatches {
    BoundaryPatches bcs;
    bcs.reserve(boundary_names.size());

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
                fmt::format("Coudln't find definition for boundary patch '{}' in boundary "
                            "conditions file '{}'",
                            bname,
                            path.string()));
        }

        auto type {doc[bname]["type"].value<std::string_view>()};

        if (!type) {
            throw std::runtime_error(fmt::format(
                "Boundary conditions for patch '{}' does not have a type in boundary condition "
                "file '{}'",
                bname,
                path.string()));
        }

        auto bc_data {parse_boundary_data(doc, bname, boundary_type_str_to_enum(type.value()))};

        bcs.emplace_back(std::string(bname), std::move(bc_data));
    }

    return bcs;
}
} // namespace prism::mesh