#pragma once

#include "pressure.h"

namespace prism::field::boundary::scalar {
// Note: This specialization is not necessary, but it is used as a demonstration of how to
// implement a custom boundary condition per field type.
template <>
class Symmetry<Pressure> : public IScalarBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(const IScalar& field, const mesh::Face& face) -> double override;
};
} // namespace prism::field::boundary::scalar
