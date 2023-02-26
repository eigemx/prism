#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <variant>
#include <vector>

#include "../types.h"

namespace prism::mesh {

enum BoundaryPatchType {
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Empty,
    Gradient,
    Unknown, // for error handling
};

struct WallBoundaryData {
    std::optional<Vector3d> velocity {std::nullopt};
    std::optional<double> temperature {std::nullopt};
};

struct InletBoundaryData {
    std::optional<Vector3d> velocity {std::nullopt};
    std::optional<double> temperature {std::nullopt};
};

struct OutletBoundaryData {
    double pressure {0.0};
};

struct GradientBoundaryData {
    std::optional<double> velocity_gradient {std::nullopt};
    std::optional<double> pressure_gradient {std::nullopt};
    std::optional<double> heat_flux {std::nullopt};
};

struct SymmetryBoundaryData {};
struct EmptyBoundaryData {};

using BoundaryData = std::variant<WallBoundaryData,
                                  InletBoundaryData,
                                  OutletBoundaryData,
                                  GradientBoundaryData,
                                  SymmetryBoundaryData,
                                  EmptyBoundaryData>;

class BoundaryPatch {
  public:
    BoundaryPatch() = delete;
    // constructor without faces
    BoundaryPatch(std::string name, BoundaryData data)
        : _name(std::move(name)), _data(std::move(data)) {
        _type = infer_boundary_type(_data);
    };

    [[nodiscard]] auto name() const -> const std::string& { return _name; }
    [[nodiscard]] auto type() const -> BoundaryPatchType { return _type; }
    [[nodiscard]] auto data() const -> const BoundaryData& { return _data; }

  private:
    auto static infer_boundary_type(const BoundaryData& data) -> BoundaryPatchType;
    std::string _name;
    BoundaryPatchType _type;
    BoundaryData _data;
};

using BoundaryPatches = std::vector<BoundaryPatch>;

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> BoundaryPatches;

} // namespace prism::mesh