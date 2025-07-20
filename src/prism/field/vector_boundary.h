#pragma once

#include "ifield.h"
#include "prism/boundary.h"
#include "prism/mesh/face.h"

namespace prism::field::boundary::vector {
class IVectorBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(const IVector& field, const mesh::Face& face) -> Vector3d = 0;
    virtual auto flux(const IVector& field, const mesh::Face& face) -> double;
    virtual auto isDirichlet() const noexcept -> bool;
};

template <IVectorBased Field>
class Fixed : public IVectorBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(const IVector& field, const mesh::Face& face) -> Vector3d override;
    auto isDirichlet() const noexcept -> bool override;
};

template <IVectorBased Field>
class NoSlip : public IVectorBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IVector& field, const mesh::Face& face) -> Vector3d override;
    auto flux(const IVector& field, const mesh::Face& face) -> double override;
    auto isDirichlet() const noexcept -> bool override;
};

template <IVectorBased Field>
class Symmetry : public IVectorBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const IVector& field, const mesh::Face& face) -> Vector3d override;
    auto flux(const IVector& field, const mesh::Face& face) -> double override;
};

template <IVectorBased Field>
class Outlet : public IVectorBoundaryHandler {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(const IVector& field, const mesh::Face& face) -> Vector3d override;
};

template <IVectorBased Field>
class ZeroGradient : public IVectorBoundaryHandler {
  public:
    auto name() const -> std::string override { return "zero-gradient"; }
    auto get(const IVector& field, const mesh::Face& face) -> Vector3d override;
};

auto inline IVectorBoundaryHandler::flux(const IVector& field, const mesh::Face& face) -> double {
    return field.valueAtFace(face).dot(face.areaVector());
}

auto inline IVectorBoundaryHandler::isDirichlet() const noexcept -> bool {
    return false;
}

template <IVectorBased Field>
auto Fixed<Field>::get(const IVector& field, const mesh::Face& face) -> Vector3d {
    const auto& patch = field.mesh()->boundaryPatch(face);
    return patch.getVectorBoundaryCondition(field.name());
}

template <IVectorBased Field>
auto Fixed<Field>::isDirichlet() const noexcept -> bool {
    return true;
}

template <IVectorBased Field>
auto NoSlip<Field>::get(const IVector& field, const mesh::Face& face) -> Vector3d {
    Fixed<Field> fixed;
    return fixed.get(field, face);
}

template <IVectorBased Field>
auto NoSlip<Field>::flux(const IVector& field, const mesh::Face& face) -> double { // NOLINT
    return 0.0; // No-slip condition implies zero flux
}

template <IVectorBased Field>
auto NoSlip<Field>::isDirichlet() const noexcept -> bool {
    return true;
}

template <IVectorBased Field>
auto Symmetry<Field>::get(const IVector& field, const mesh::Face& face) -> Vector3d {
    Vector3d Uc = field.valueAtCell(face.owner());
    const Vector3d& normal = face.normal();

    // Project out the normal component to enforce zero normal velocity at symmetry plane
    // The projection of a vector 'v' onto a normal 'n' is (v.dot(n)) * n
    // So, the tangential component is v - (v.dot(n)) * n
    return Uc - Uc.dot(normal) * normal;
}

template <IVectorBased Field>
auto Symmetry<Field>::flux(const IVector& field, const mesh::Face& face) -> double { // NOLINT
    return 0.0;
}

template <IVectorBased Field>
auto Outlet<Field>::get(const IVector& field, const mesh::Face& face) -> Vector3d {
    ZeroGradient<Field> zg;
    return zg.get(field, face);
}

template <IVectorBased Field>
auto ZeroGradient<Field>::get(const IVector& field, const mesh::Face& face) -> Vector3d {
    return field.valueAtCell(face.owner());
}

} // namespace prism::field::boundary::vector
