#include "velocity.h"

namespace prism::field {
VelocityComponent::VelocityComponent(std::string name,
                                     const mesh::PMesh& mesh,
                                     double value,
                                     Coord coord,
                                     IVector* parent)
    : Scalar(std::move(name), mesh, value, parent), _coord(coord) {}

VelocityComponent::VelocityComponent(std::string name,
                                     const mesh::PMesh& mesh,
                                     VectorXd data,
                                     Coord coord,
                                     IVector* parent)
    : Scalar(std::move(name), mesh, std::move(data), parent), _coord(coord) {}

VelocityComponent::VelocityComponent(std::string name,
                                     const mesh::PMesh& mesh,
                                     VectorXd data,
                                     VectorXd face_data,
                                     Coord coord,
                                     IVector* parent)
    : Scalar(std::move(name), mesh, std::move(data), std::move(face_data), parent),
      _coord(coord) {}

Velocity::Velocity(std::string name, const mesh::PMesh& mesh, double value)
    : Vector(std::move(name), mesh, value) {}

Velocity::Velocity(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : Vector(std::move(name), mesh, data) {}

Velocity::Velocity(std::string name,
                   const mesh::PMesh& mesh,
                   std::array<VelocityComponent, 3>& fields)
    : Vector(std::move(name), mesh, fields) {}

} // namespace prism::field