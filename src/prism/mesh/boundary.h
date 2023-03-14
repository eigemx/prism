#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <variant>
#include <vector>

#include "../types.h"

namespace prism::mesh {

enum class BoundaryPatchType {
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Empty,
    Gradient,
    Unknown, // for error handling
};

enum class BoundaryAttributeType { Scalar, Vector };

using BoundaryAttributeData = std::variant<double, Vector3d>;

class BoundaryPatchAttribute {
  public:
    BoundaryPatchAttribute(std::string name,
                           BoundaryAttributeType type,
                           BoundaryAttributeData data)
        : _name(std::move(name)), _type(type), _data(std::move(data)) {}

    [[nodiscard]] auto name() const noexcept -> const std::string& { return _name; }
    [[nodiscard]] auto type() const noexcept -> BoundaryAttributeType { return _type; }
    [[nodiscard]] auto data() const noexcept -> const BoundaryAttributeData& { return _data; }

  private:
    std::string _name;
    BoundaryAttributeType _type;
    BoundaryAttributeData _data;
};

using BoundaryPatchAttributes = std::vector<BoundaryPatchAttribute>;

class BoundaryPatch {
  public:
    BoundaryPatch() = delete;

    BoundaryPatch(std::string name, BoundaryPatchAttributes attributes, BoundaryPatchType type)
        : _name(std::move(name)), _attributes(std::move(attributes)), _type(type) {};

    [[nodiscard]] auto name() const noexcept -> const std::string& { return _name; }
    [[nodiscard]] auto type() const noexcept -> BoundaryPatchType { return _type; }
    [[nodiscard]] auto attributes() const noexcept -> const BoundaryPatchAttributes& {
        return _attributes;
    }

  private:
    std::string _name;
    BoundaryPatchType _type;
    BoundaryPatchAttributes _attributes;
};

using BoundaryPatches = std::vector<BoundaryPatch>;

auto read_boundary_conditions(const std::filesystem::path& path,
                              const std::vector<std::string_view>& boundary_names)
    -> BoundaryPatches;

} // namespace prism::mesh