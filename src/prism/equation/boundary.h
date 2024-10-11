#pragma once

#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/field/velocity.h"
#include "prism/mesh/boundary.h"

namespace prism::eqn {

// forward declaration
template <field::IScalarBased Field>
class Transport;

using Momentum = Transport<field::VelocityComponent>;

} // namespace prism::eqn

namespace prism::eqn::boundary {

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

template <>
class NoSlip<Momentum> : public IEquationBoundaryHandler<Momentum> {
  public:
    auto name() const noexcept -> std::string override { return "no-slip"; }
    void apply(Momentum& eqn, const mesh::BoundaryPatch& patch) override;
};

} // namespace prism::eqn::boundary
