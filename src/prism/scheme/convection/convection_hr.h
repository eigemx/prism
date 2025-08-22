#pragma once

#include "prism/field/ifield.h"
#include "prism/mesh/cell.h"

namespace prism::scheme::convection {

auto phiAtDummyUpwind(SharedPtr<field::IScalar>& field,
                      const mesh::Cell& cell,
                      const mesh::Cell& downwind) -> f64;

auto phiTilde(SharedPtr<field::IScalar> field, const mesh::Cell& cell, const mesh::Cell& downwind)
    -> f64;

} // namespace prism::scheme::convection
