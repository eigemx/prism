#include <string>

#include "prism/boundary.h"
#include "prism/field/field.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::gradient::boundary {
class GradSchemeBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d = 0;
};

class Empty : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "empty"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};

class Symmetry : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};

class Outlet : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};

class Fixed : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};

class VelocityInlet : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};

class FixedGradient : public GradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed=gradient"; }
    auto get(const prism::field::Scalar& field, const prism::mesh::Face& face)
        -> prism::Vector3d override;
};


} // namespace prism::gradient::boundary