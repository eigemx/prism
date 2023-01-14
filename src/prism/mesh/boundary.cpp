#include "boundary.h"

#include <fmt/format.h>
#include <toml++/toml.h>

#include <fstream>
#include <iostream>
#include <string_view>


namespace prism::mesh {
auto boundary_type_str_to_enum(std::string_view type) -> BoundaryConditionType {
    if (type == "wall") {
        return BoundaryConditionType::Wall;
    }
    if (type == "inlet") {
        return BoundaryConditionType::Inlet;
    }
    if (type == "outlet") {
        return BoundaryConditionType::Outlet;
    }
    if (type == "gradient") {
        return BoundaryConditionType::Gradient;
    }
    if (type == "symmetry") {
        return BoundaryConditionType::Symmetry;
    }
    if (type == "empty") {
        return BoundaryConditionType::Empty;
    }

    return BoundaryConditionType::Unknown;
}

auto parse_wall(const toml::table& doc, std::string_view bname) -> WallBoundaryData {
    auto velocity = doc[bname]["velocity"];
    auto temperature = doc[bname]["temperature"];

    if (!velocity && !temperature) {
        throw std::runtime_error(fmt::format(
            "Boundary '{}' is defined as a wall and has no velocity or temperature", bname));
    }

    std::optional<Vector3d> velocity_vec {std::nullopt};
    std::optional<double> temperature_val {std::nullopt};

    if (velocity) {
        // check if velocity is a vector
        if (!velocity.is_array() || velocity.as_array()->size() != 3) {
            throw std::runtime_error(
                fmt::format("Wall boundary '{}' has a velocity that is not a 3D vector", bname));
        }

        // check if velocity is a vector of numbers
        for (std::size_t i = 0; i < 3; ++i) {
            if (!velocity[i].is_floating_point() || !velocity[i].is_number()) {
                throw std::runtime_error(fmt::format(
                    "Wall boundary '{}' has a velocity that is not a 3D vector of numbers",
                    bname));
            }
        }

        velocity_vec =
            Vector3d {velocity[0].value<double>().value(), velocity[1].value<double>().value(),
                      velocity[2].value<double>().value()};
    }

    if (temperature) {
        if (!temperature.is_floating_point() || !temperature.is_number()) {
            throw std::runtime_error(
                fmt::format("Wall boundary '{}' has a temperature that is not a number", bname));
        }

        temperature_val = temperature.value<double>().value();
    }

    return WallBoundaryData {velocity_vec, temperature_val};
}

auto parse_inlet(const toml::table& doc, std::string_view bname) -> InletBoundaryData {
    auto velocity = doc[bname]["velocity"];
    auto temperature = doc[bname]["temperature"];

    if (!velocity && !temperature) {
        throw std::runtime_error(fmt::format(
            "Boundary '{}' is defined as an inlet and has no velocity or temperature", bname));
    }

    std::optional<Vector3d> velocity_vec {std::nullopt};
    std::optional<double> temperature_val {std::nullopt};

    if (velocity) {
        // check if velocity is a vector
        if (!velocity.is_array() || velocity.as_array()->size() != 3) {
            throw std::runtime_error(
                fmt::format("Inlet boundary '{}' has a velocity that is not a 3D vector", bname));
        }

        // check if velocity is a vector of numbers
        for (std::size_t i = 0; i < 3; ++i) {
            if (!velocity[i].is_floating_point()) {
                throw std::runtime_error(fmt::format(
                    "Inlet boundary '{}' has a velocity that is not a 3D vector of numbers",
                    bname));
            }
        }

        velocity_vec =
            Vector3d {velocity[0].value<double>().value(), velocity[1].value<double>().value(),
                      velocity[2].value<double>().value()};
    }

    if (temperature) {
        if (!temperature.is_floating_point() || !temperature.is_number()) {
            throw std::runtime_error(
                fmt::format("Inlet boundary '{}' has a temperature that is not a number", bname));
        }

        temperature_val = temperature.value<double>().value();
    }

    return InletBoundaryData {velocity_vec, temperature_val};
}

auto parse_outlet(const toml::table& doc, std::string_view bname) -> OutletBoundaryData {
    auto pressure = doc[bname]["pressure"];

    if (!pressure) {
        throw std::runtime_error(
            fmt::format("Boundary '{}' is defined as an outlet and has no pressure", bname));
    }

    if (!pressure.is_floating_point()) {
        throw std::runtime_error(
            fmt::format("Outlet boundary '{}' has a pressure that is not a number", bname));
    }

    return OutletBoundaryData {pressure.value<double>().value()};
}

auto parse_gradient(const toml::table& doc, std::string_view bname) -> GradientBoundaryData {
    // check velocity_gradient, pressure_gradient, and heat flux
    auto velocity_gradient = doc[bname]["velocity_gradient"];
    auto pressure_gradient = doc[bname]["pressure_gradient"];
    auto heat_flux = doc[bname]["heat_flux"];

    if (!velocity_gradient && !pressure_gradient && !heat_flux) {
        throw std::runtime_error(fmt::format(
            "Boundary '{}' is defined as a gradient and has no velocity gradient, pressure "
            "gradient, or heat flux",
            bname));
    }

    std::optional<double> velocity_gradient_val {std::nullopt};
    std::optional<double> pressure_gradient_val {std::nullopt};
    std::optional<double> heat_flux_val {std::nullopt};

    if (velocity_gradient) {
        if (!velocity_gradient.is_floating_point() || !velocity_gradient.is_number()) {
            throw std::runtime_error(fmt::format(
                "Gradient boundary '{}' has a velocity gradient that is not a number", bname));
        }

        velocity_gradient_val = velocity_gradient.value<double>().value();
    }

    if (pressure_gradient) {
        if (!pressure_gradient.is_floating_point() || !pressure_gradient.is_number()) {
            throw std::runtime_error(fmt::format(
                "Gradient boundary '{}' has a pressure gradient that is not a number", bname));
        }

        pressure_gradient_val = pressure_gradient.value<double>().value();
    }

    if (heat_flux) {
        if (!heat_flux.is_floating_point() || !heat_flux.is_number()) {
            throw std::runtime_error(fmt::format(
                "Gradient boundary '{}' has a heat flux that is not a number", bname));
        }

        heat_flux_val = heat_flux.value<double>().value();
    }

    return GradientBoundaryData {velocity_gradient_val, pressure_gradient_val, heat_flux_val};
}

auto parse_boundary_data(const toml::table& table, std::string_view bname,
                         BoundaryConditionType bc_type) -> BoundaryData {
    switch (bc_type) {
        case BoundaryConditionType::Wall:
            return parse_wall(table, bname);

        case BoundaryConditionType::Inlet:
            return parse_inlet(table, bname);

        case BoundaryConditionType::Outlet:
            return parse_outlet(table, bname);

        case BoundaryConditionType::Gradient:
            return parse_gradient(table, bname);

        case BoundaryConditionType::Symmetry:
            return SymmetryBoundaryData {};

        case BoundaryConditionType::Empty:
            return EmptyBoundaryData {};

        case BoundaryConditionType::Unknown:
            throw std::runtime_error("Unknown boundary condition type");
    }
}

auto BoundaryCondition::infer_boundary_type(const BoundaryData& data) -> BoundaryConditionType {
    if (std::holds_alternative<WallBoundaryData>(data)) {
        return BoundaryConditionType::Wall;
    }
    if (std::holds_alternative<InletBoundaryData>(data)) {
        return BoundaryConditionType::Inlet;
    }
    if (std::holds_alternative<OutletBoundaryData>(data)) {
        return BoundaryConditionType::Outlet;
    }
    if (std::holds_alternative<GradientBoundaryData>(data)) {
        return BoundaryConditionType::Gradient;
    }
    if (std::holds_alternative<SymmetryBoundaryData>(data)) {
        return BoundaryConditionType::Symmetry;
    }
    if (std::holds_alternative<EmptyBoundaryData>(data)) {
        return BoundaryConditionType::Empty;
    }
    throw std::runtime_error("Non-implemented boundary type");
}

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> BoundaryConditions {
    BoundaryConditions bcs;
    bcs.reserve(boundary_names.size());

    auto fstream {std::ifstream {path}};
    toml::table doc;

    try {
        doc = toml::parse(fstream);
    } catch (const toml::parse_error& e) {
        throw std::runtime_error(
            fmt::format("Failed to parse boundary condition file: '{}', complete error: {}",
                        path.string(), e.what()));
    }

    for (const auto& bname : boundary_names) {
        auto table {doc[bname.data()]};
        if (!table) {
            throw std::runtime_error(
                fmt::format("Boundary patch '{}' not found in boundary condition file {}", bname,
                            path.string()));
        }

        auto type {doc[bname]["type"].value<std::string_view>()};

        if (!type) {
            throw std::runtime_error(fmt::format(
                "Boundary conditions for patch '{}' does not have a type in boundary condition "
                "file {}",
                bname, path.string()));
        }

        auto bc_data {parse_boundary_data(doc, bname, boundary_type_str_to_enum(type.value()))};

        bcs.emplace_back(std::string(bname), std::move(bc_data));
    }

    return bcs;
}
} // namespace prism::mesh