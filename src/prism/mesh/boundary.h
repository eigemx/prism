#pragma once

#include <string>
#include <variant>
#include <vector>

#include "../types.h"

namespace prism::mesh {

enum BoundaryType {
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Empty,
    Gradient,
};

struct WallBoundaryData {
    Vector3d velocity;
    double temperature;
};

struct InletBoundaryData {
    Vector3d velocity;
    double temperature;
};

struct OutletBoundaryData {
    double pressure;
};

struct GradientBoundaryData {
    double velocity_gradient;
    double pressure_gradient;
    double temperature_gradient;
};

using BoundaryData =
    std::variant<WallBoundaryData, InletBoundaryData, OutletBoundaryData, GradientBoundaryData>;

class BoundaryPatch {
public:
    BoundaryPatch() = delete;

private:
    std::string _name;
    std::vector<std::size_t> _faces;
    BoundaryType type;
    BoundaryData data;
};


} // namespace prism::mesh