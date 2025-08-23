#pragma once

#include "prism/field/scalar.h"
#include "prism/mesh/cell.h"

namespace prism::scheme::convection {

auto phiAtDummyUpwind(SharedPtr<field::Scalar>& field,
                      const mesh::Cell& cell,
                      const mesh::Cell& downwind) -> f64;

auto phiTilde(SharedPtr<field::Scalar>& field, const mesh::Cell& cell, const mesh::Cell& downwind)
    -> f64;

} // namespace prism::scheme::convection
