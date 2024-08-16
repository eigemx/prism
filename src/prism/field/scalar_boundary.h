#pragma once

#include "boundary.h"
#include "prism/mesh/face.h"
#include "scalar.h"


namespace prism::field::boundary {

template <>
class Fixed<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class NoSlip<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class VelocityInlet<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class Empty<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "empty"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class Symmetry<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class Outlet<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

template <>
class FixedGradient<Scalar> : public FieldBoundaryHandler<Scalar> {
  public:
    auto name() const -> std::string override { return "fixed-gradient"; }
    auto get(const Scalar& field, const mesh::Face& face) -> double override;
};

} // namespace prism::field::boundary