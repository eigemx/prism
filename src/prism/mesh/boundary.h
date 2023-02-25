#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "../types.h"

namespace prism::mesh {

enum BoundaryConditionType {
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

class BoundaryCondition {
  public:
    BoundaryCondition() = delete;
    // constructor without faces
    BoundaryCondition(std::string&& name, BoundaryData&& data)
        : _name(std::move(name)), _data(std::move(data)) {
        _type = infer_boundary_type(_data);
    };

    BoundaryCondition(std::string&& name, std::vector<std::size_t>&& faces, BoundaryData&& data)
        : _name(std::move(name)), _faces(std::move(faces)), _data(std::move(data)) {
        _type = infer_boundary_type(_data);
    };

    [[nodiscard]] auto name() const -> const std::string& { return _name; }
    [[nodiscard]] auto faces() const -> const std::vector<std::size_t>& { return _faces; }
    [[nodiscard]] auto type() const -> BoundaryConditionType { return _type; }
    [[nodiscard]] auto data() const -> const BoundaryData& { return _data; }

    void set_faces(std::vector<std::size_t>&& faces) { _faces = std::move(faces); }


  private:
    auto static infer_boundary_type(const BoundaryData& data) -> BoundaryConditionType;
    std::string _name;
    std::vector<std::size_t> _faces;
    BoundaryConditionType _type;
    BoundaryData _data;
};

using BoundaryConditions = std::vector<BoundaryCondition>;

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> BoundaryConditions;

} // namespace prism::mesh