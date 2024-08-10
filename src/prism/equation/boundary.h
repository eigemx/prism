#pragma once

#include <prism/field/velocity.h>

#include "prism/boundary.h"
#include "prism/mesh/boundary.h"

namespace prism {

// forward declaration
template <typename Field>
class TransportEquation;

using MomentumEquation = TransportEquation<field::VelocityComponent>;

} // namespace prism

namespace prism::boundary {

template <typename Equation>
class IEquationBoundaryHandler : public IBoundaryHandler {
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
class NoSlip<MomentumEquation> : public IEquationBoundaryHandler<MomentumEquation> {
  public:
    auto name() const noexcept -> std::string override { return "no-slip"; }
    void apply(MomentumEquation& eqn, const mesh::BoundaryPatch& patch) override;
};


} // namespace prism::boundary