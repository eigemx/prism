#pragma once

#include "ifield.h"
#include "prism/boundary.h"
#include "prism/mesh/face.h"

namespace prism::field::boundary {
class IScalarBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const IScalar& field, const mesh::Face& face) -> double = 0;
};

template <IScalarBased Field>
class Fixed : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class NoSlip : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class VelocityInlet : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class Empty : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "empty"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class Symmetry : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class Outlet : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class FixedGradient : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed-gradient"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
class ZeroGradient : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "zero-gradient"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
auto Fixed<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh().boundaryPatch(face);
    return patch.getScalarBoundaryCondition(field.name());
}

template <IScalarBased Field>
auto Empty<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    // TODO: Empty faces field value should not contribute to the solution, we need to avoid
    // having an "Empty" handler for fields.
    return field.valueAtCell(face.owner());
}

template <IScalarBased Field>
auto Symmetry<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    ZeroGradient<Field> zg;
    return zg.get(field, face);
}

template <IScalarBased Field>
auto Outlet<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    ZeroGradient<Field> zg;
    return zg.get(field, face);
}

template <IScalarBased Field>
auto FixedGradient<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    const auto& mesh = field.mesh();
    const auto& patch = mesh.boundaryPatch(face);

    Vector3d grad_at_boundary = patch.getVectorBoundaryCondition(field.name());
    const auto& owner = mesh.cell(face.owner());

    Vector3d e = face.center() - owner.center();
    double d_Cf = e.norm();
    e = e / e.norm();
    grad_at_boundary = grad_at_boundary * d_Cf;

    return grad_at_boundary.dot(e) + field.valueAtCell(owner);
}

template <IScalarBased Field>
auto ZeroGradient<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    return field.valueAtCell(face.owner());
}

} // namespace prism::field::boundary