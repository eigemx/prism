#pragma once

#include "ifield.h"
#include "prism/boundary.h"
#include "prism/mesh/face.h"

namespace prism::field::boundary::scalar {
class IScalarBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const IScalar& field, const mesh::Face& face) -> double = 0;
    virtual auto isDirichlet() const noexcept -> bool { return false; }
};

template <IScalarBased Field>
class Fixed : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
    auto isDirichlet() const noexcept -> bool override { return true; }
};

template <IScalarBased Field>
class NoSlip : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
    auto isDirichlet() const noexcept -> bool override { return true; }
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
class ZeroGradient : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "zero-gradient"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <IScalarBased Field>
auto Fixed<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    const auto& patch = field.mesh()->boundaryPatch(face);
    return patch.getScalarBoundaryCondition(field.name());
}

template <IScalarBased Field>
auto NoSlip<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    Fixed<Field> fixed;
    return fixed.get(field, face);
}

template <IScalarBased Field>
auto Symmetry<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    return field.valueAtCell(face.owner());
}

template <IScalarBased Field>
auto Outlet<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    return field.valueAtCell(face.owner());
}

template <IScalarBased Field>
auto ZeroGradient<Field>::get(const IScalar& field, const mesh::Face& face) -> double {
    /// TODO: is this correct? this is based on 15.160, but I am note sure if it's valid for zero
    /// gradient boundary condition.
    const auto& owner = field.mesh()->cell(face.owner());
    const Vector3d d_Cb = face.center() - owner.center();
    return field.valueAtCell(owner); // + field.gradAtCellStored(owner).dot(d_Cb);
}

} // namespace prism::field::boundary::scalar
