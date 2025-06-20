#pragma once

#include "pressure.h"

namespace prism::field::boundary {
/// TODO: remove this as we don't need a NoSlip boundary condition for pressure fields. NoSlip is
/// applicable for velocity fields, where pressure in such boundaries is set as zero gradient.
template <>
class NoSlip<Pressure> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};
} // namespace prism::field::boundary