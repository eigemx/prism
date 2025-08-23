#pragma once

#include "prism/boundary.h"
#include "prism/field/velocity.h"
#include "prism/mesh/boundary.h"
#include "prism/scheme/scheme.h"

namespace prism::eqn {

// forward declaration
class Transport;

class Momentum;

} // namespace prism::eqn

namespace prism::eqn::boundary {

template <typename To>
auto castScheme(const SharedPtr<prism::scheme::IScheme>& ptr) -> SharedPtr<To> {
    return std::dynamic_pointer_cast<To>(ptr);
}

template <typename Equation>
class IEquationBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const noexcept -> std::string = 0;
    virtual void apply(Equation& eqn, const mesh::BoundaryPatch& patch) = 0;
};

template <typename Equation>
class NoSlip : public IEquationBoundaryHandler<Equation> {
  public:
    auto name() const noexcept -> std::string override { return "no-slip"; }
    void apply(Equation& eqn, const mesh::BoundaryPatch& patch) override;
};

template <typename Equation>
class Symmetry : public IEquationBoundaryHandler<Equation> {
  public:
    auto name() const noexcept -> std::string override { return "symmetry"; }
    void apply(Equation& eqn, const mesh::BoundaryPatch& patch) override;
};

template <typename Equation>
class Outlet : public IEquationBoundaryHandler<Equation> {
  public:
    auto name() const noexcept -> std::string override { return "outlet"; }
    void apply(Equation& eqn, const mesh::BoundaryPatch& patch) override;
};

template <>
class NoSlip<Momentum> : public IEquationBoundaryHandler<Momentum> {
  public:
    auto name() const noexcept -> std::string override { return "no-slip"; }
    void apply(Momentum& eqn, const mesh::BoundaryPatch& patch) override;
};

template <>
class Symmetry<Momentum> : public IEquationBoundaryHandler<Momentum> {
  public:
    auto name() const noexcept -> std::string override { return "symmetry"; }
    void apply(Momentum& eqn, const mesh::BoundaryPatch& patch) override;
};

template <>
class Outlet<Momentum> : public IEquationBoundaryHandler<Momentum> {
  public:
    auto name() const noexcept -> std::string override { return "outlet"; }
    void apply(Momentum& eqn, const mesh::BoundaryPatch& patch) override;
};

} // namespace prism::eqn::boundary
