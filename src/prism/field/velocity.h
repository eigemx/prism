#pragma once

#include <array>
#include <optional>
#include <string>

#include "scalar.h"
#include "units.h"
#include "vector.h"

namespace prism::field {
class Velocity;

class VelocityComponent : public Scalar, public units::VelocityUnit {
  public:
    VelocityComponent(std::string name,
                      const mesh::PMesh& mesh,
                      double value,
                      Coord coord,
                      IVector* parent = nullptr);
    VelocityComponent(std::string name,
                      const mesh::PMesh& mesh,
                      VectorXd data,
                      Coord coord,
                      IVector* parent = nullptr);
    VelocityComponent(std::string name,
                      const mesh::PMesh& mesh,
                      VectorXd data,
                      VectorXd face_data,
                      Coord coord,
                      IVector* parent = nullptr);

    auto coord() const noexcept -> std::optional<Coord> override { return _coord; }

  private:
    Coord _coord;
};

class Velocity : public detail::Vector<VelocityComponent>, public units::VelocityUnit {
  public:
    Velocity(std::string name, const mesh::PMesh& mesh, double value);
    Velocity(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    Velocity(std::string name, const mesh::PMesh& mesh, std::array<VelocityComponent, 3>& fields);
};

} // namespace prism::field