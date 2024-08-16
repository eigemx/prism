#pragma once

#include "pressure.h"

namespace prism::field::boundary {
template <>
class NoSlip<Pressure> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};
} // namespace prism::field::boundary