#pragma once

#include "prism/boundary.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/pmesh.h"

namespace prism::field {
class Scalar;
}

// In convection & diffusion schemes, we defined the boundary condition handlers directly in
// namespace prism::boundary, but since field.h is included in most header files there will be
// redefinition of type conflict (when defining boundary handlers again for diffusion &
// convection), so we avoid this by introducing a nested namespace inside prism::field to include
// the definition of field boundary handlers.
namespace prism::field::boundary {

template <typename Field, typename ValueType>
class FieldBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const field::Scalar& field, const mesh::Face& face) -> ValueType = 0;
};

class Fixed : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

class VelocityInlet : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

class Empty : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "empty"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

class Symmetry : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

class Outlet : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

class FixedGradient : public FieldBoundaryHandler<field::Scalar, double> {
  public:
    auto name() const -> std::string override { return "fixed-gradient"; }
    auto get(const field::Scalar& field, const mesh::Face& face) -> double override;
};

} // namespace prism::field::boundary