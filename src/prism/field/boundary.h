#pragma once

#include "prism/boundary.h"
#include "prism/mesh/face.h"


namespace prism::field::boundary {

template <typename Field>
class FieldBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const Field& field, const mesh::Face& face) -> double = 0;
};

template <typename Field>
class Fixed : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};


template <typename Field>
class NoSlip : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

template <typename Field>
class VelocityInlet : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

template <typename Field>
class Empty : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "empty"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

template <typename Field>
class Symmetry : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

template <typename Field>
class Outlet : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

template <typename Field>
class FixedGradient : public FieldBoundaryHandler<Field> {
  public:
    auto name() const -> std::string override { return "fixed-gradient"; }
    auto get(const Field& field, const mesh::Face& face) -> double override;
};

} // namespace prism::field::boundary