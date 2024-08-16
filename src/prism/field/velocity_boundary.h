#pragma once

#include "boundary.h"
#include "velocity.h"

namespace prism::field::boundary {
template <>
class NoSlip<VelocityComponent> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};

template <>
class VelocityInlet<VelocityComponent> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "velocity-inlet"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};
} // namespace prism::field::boundary