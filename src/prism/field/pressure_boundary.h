#pragma once

#include "pressure.h"

namespace prism::field::boundary::scalar {
template <>
class Symmetry<Pressure> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};
} // namespace prism::field::boundary::scalar